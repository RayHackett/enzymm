import argparse
import collections
from pathlib import Path
from typing import List, Tuple, Dict, Optional, IO
from functools import cached_property
from dataclasses import dataclass, field
import warnings
import csv
import itertools
import io
import sys
import os

import pyjess

from .template import AnnotatedTemplate, Vec3, load_templates, check_template
from .utils import chunks, ranked_argsort

# Matcher goes Brrrrr
from multiprocessing.pool import ThreadPool
import functools

__all__ = [
    "Match",
    "single_query_run",
    "Matcher",
]

@dataclass
class Match:
    """"Class for storing annotated Jess hits. Wrapper around pyjess.Hit and the original Template object"""
    hit: pyjess.Hit
    template: AnnotatedTemplate
    complete: bool = field(default=False)
    index: int = field(default=0)

    def dumps(self, header:bool = False) -> str:
        buffer = io.StringIO()
        self.dump(buffer, header=header)
        return buffer.getvalue() # returns entire content temporary file object as a string

    def dump2pdb(self, file: IO[str], include_query:bool = False, transform:bool = False):

        def write_atom_line(atom: pyjess.Atom) -> str:
            one_char_elements = {'H', 'B', 'C', 'N', 'O', 'F', 'P', 'S', 'K', 'V', 'Y', 'I', 'W', 'U'}
            if atom.element in one_char_elements:
                return f"ATOM  {atom.serial:>5}  {atom.name:<3s}{atom.altloc if atom.altloc is not None else '':<1}{atom.residue_name:<3}{atom.chain_id:>2}{atom.residue_number:>4}{atom.insertion_code:1s}   {atom.x:>8.3f}{atom.y:>8.3f}{atom.z:>8.3f}{atom.occupancy:>6.2f}{atom.temperature_factor:>6.2f}      {atom.segment:<4s}{atom.element:>2s} \n"
            else:
                return f"ATOM  {atom.serial:>5} {atom.name:<4s}{atom.altloc if atom.altloc is not None else '':<1}{atom.residue_name:<3}{atom.chain_id:>2}{atom.residue_number:>4}{atom.insertion_code:1s}   {atom.x:>8.3f}{atom.y:>8.3f}{atom.z:>8.3f}{atom.occupancy:>6.2f}{atom.temperature_factor:>6.2f}      {atom.segment:<4s}{atom.element:>2s} \n"

        if include_query: # write the original query molecule too
            file.write(f'REMARK MOLECULE_ID {self.hit.molecule.id}\n')
            for atom in self.hit.molecule:
                file.write(write_atom_line(atom))
            file.write('END\n\n')

        file.write(f'REMARK TEMPLATE_PDB {str(self.template.pdb_id)}_{",".join(set(res.chain_id for res in self.template.residues))}\n')
        if self.template.cluster:
            file.write(f'REMARK TEMPLATE CLUSTER {str(self.template.cluster.id)}_{str(self.template.cluster.member)}_{str(self.template.cluster.size)}\n')
        if self.template.represented_sites:
            file.write(f'REMARK TEMPLATE RESIDUES {str(self.template.template_id_string)}\n')
        file.write(f'REMARK MOLECULE_ID {str(self.hit.molecule.id)}\n')
        file.write(f'MATCH INDEX {self.index}\n')

        if transform:
            file.write('REMARK TEMPLATE COORDINATE FRAME\n')
        else:
            file.write('REMARK QUERY COORDINATE FRAME\n')

        for atom in self.hit.atoms(transform=transform): # if transform == True then the atom coordinates are transformed to the template reference frame
            file.write(write_atom_line(atom))
        file.write('END\n\n')

    def dump(self, file: IO[str], header: bool = False):
        writer = csv.writer(file, dialect="excel-tab", delimiter='\t', lineterminator='\n')
        if header:
                writer.writerow([
                    "query_id",
                    "match_index",
                    "template_pdb_id",
                    "template_pdb_chains",
                    "template_cluster_id",
                    "template_cluster_member",
                    "template_cluster_size",
                    "template_effective_size",
                    "template_dimension",
                    "template_mcsa_id",
                    "template_uniprot_id",
                    "template_ec",
                    "template_cath",
                    "template_multimeric",
                    "query_multimeric",
                    "query_atom_count",
                    "query_residue_count",
                    "rmsd",
                    "log_evalue",
                    "orientation",
                    "preserved_order",
                    "completeness",
                    "matched_residues"])

        writer.writerow([
                str(self.hit.molecule.id),
                str(self.index),
                str(self.template.pdb_id if self.template.pdb_id else ''),
                (','.join(set(res.chain_id for res in self.template.residues))),
                str(self.template.cluster.id if self.template.cluster else ''),
                str(self.template.cluster.member if self.template.cluster else ''),
                str(self.template.cluster.size if self.template.cluster else ''),
                str(self.template.effective_size),
                str(self.template.dimension),
                str(self.template.mcsa_id if self.template.mcsa_id else ''),
                str(self.template.uniprot_id if self.template.uniprot_id else ''),
                ",".join(self.template.ec if self.template.ec is not None else ''),
                ",".join(self.template.cath if self.template.cath else ''),
                str(self.template.multimeric),
                str(self.multimeric),
                str(self.query_atom_count),
                str(self.query_residue_count),
                str(self.hit.rmsd),
                str(self.hit.log_evalue),
                str(self.orientation),
                str(self.preserved_resid_order),
                str(self.complete),
                (','.join('_'.join(t) for t in self.matched_residues))])

    @cached_property
    def atom_triplets(self) -> List[Tuple[pyjess.Atom, pyjess.Atom, pyjess.Atom]]:
        # list with matched residues
        # # Hit.atoms is a list of matched atoms with all info on residue numbers and residue chain ids and atom types, this should conserve order if Hit.atoms is a list!!!
        atom_triplets = []
        for atom_triplet in chunks(self.hit.atoms(transform=True), 3): # yield chunks of 3 atoms each, transform true because for angle calculation atoms need to be in template reference frame
            if len(atom_triplet) != 3:
                raise ValueError(f'Failed to construct residues. Got only {len(atom_triplet)} ATOM lines')
            # check if all three atoms belong to the same residue by adding a tuple of their residue defining properties to a set
            unique_residues = { (atom.residue_name, atom.chain_id, atom.residue_number) for atom in atom_triplet }
            if len(unique_residues) != 1:
                raise ValueError(f'Mixed up atom triplets {unique_residues}. The atoms come from different residues!')
            atom_triplets.append(atom_triplet)
        return atom_triplets

    @property
    def matched_residues(self) -> List[Tuple[str, str, str]]:
        return [ (atom_triplet[0].residue_name, atom_triplet[0].chain_id, str(atom_triplet[0].residue_number)) for atom_triplet in self.atom_triplets ]

    @property
    def multimeric(self) -> bool:
        """`bool`: Boolean if the atoms in the hit stem from multiple protein chains
        """
        # note that these are pyjess atom objects!
        return not all(atom.chain_id == self.hit.atoms()[0].chain_id for atom in self.hit.atoms())

    @property
    def preserved_resid_order(self) -> bool:
        """`bool`: Boolean if the residues in the template and in the matched query structure have the same relative order.
        This is a good filtering parameter but excludes hits on examples of convergent evolution or circular permutations
        """
        if self.template.multimeric or self.multimeric:
            return False
        else:
            # Now extract relative atom order in hit
            return ranked_argsort([atom_triplet[0].residue_number for atom_triplet in self.atom_triplets]) == self.template.relative_order

    @cached_property
    def match_vector_list(cls) -> List[Vec3]:
        # !!! atom coordinates must be in template coordinate system!
        vector_list = []
        for residue_index, residue in enumerate(cls.template.residues):
            first_atom_index, second_atom_index = residue.orientation_vector_indices
            if second_atom_index == 9: # Calculate orientation vector going from middle_atom to mitpoint between side1 and side2
                middle_atom = cls.atom_triplets[residue_index][first_atom_index]
                side1, side2 = [atom for atom in cls.atom_triplets[residue_index] if atom != middle_atom]
                midpoint = (Vec3.from_xyz(side1) + Vec3.from_xyz(side2)) / 2
                vector_list.append(midpoint - Vec3.from_xyz(middle_atom))
            else:
                # Calculate orientation vector going from first_atom to second_atom_index
                first_atom = cls.atom_triplets[residue_index][first_atom_index]
                second_atom = cls.atom_triplets[residue_index][second_atom_index]
                vector_list.append(Vec3.from_xyz(second_atom) - Vec3.from_xyz(first_atom))
        return vector_list

    @property
    def template_vector_list(self) -> List[Vec3]:
        return [res.orientation_vector for res in self.template.residues]
    
    @property
    def orientation(self) -> float: # average angle
        if len(self.template_vector_list) != len(self.match_vector_list):
            raise ValueError('Vector lists for Template and matching Query structure had different lengths.')

        # now calculate the angle between the vector of the template and the query per residue
        angle_list = []
        for i in range(len(self.template_vector_list)):
            angle_list.append(self.template_vector_list[i].angle_to(self.match_vector_list[i]))
        return sum(angle_list)/len(angle_list)

    @property
    def query_atom_count(self) -> int:
        return len(self.hit.molecule)

    @property
    def query_residue_count(self) -> int:
        all_residue_numbers = set()
        for atom in self.hit.molecule:
            all_residue_numbers.add(atom.residue_number)
        return len(all_residue_numbers)

def single_query_run(molecule: pyjess.Molecule, templates: List[AnnotatedTemplate], rmsd_threshold: float = 2.0, distance_cutoff: float = 1.5, max_dynamic_distance: float = 1.5, max_candidates: int = 10000) -> List[Match]:
    """`list` of `Match`: Match the list of Templates to one Molecule (pyjess.Molecule) and caluclate completeness for each match (True if all members of a Template cluster match)"""

    # mapping pyjess.Template to AnnotatedTemplate via the str AnnotatedTemplate.id which is preserved in pyjess.Template.id
    id_to_template: Dict[str, AnnotatedTemplate] = {}
    for template in templates:
        if template.id in id_to_template:
            raise KeyError('Multiple Templates share the same id! To be able to map annotations, create and pass AnnoateTemplate objects with unique ids!')
        id_to_template[template.id] = template

    # killswitch is controlled by max_candidates. Internal default is currently 1000
    # killswitch serves to limit the iterations in cases where the template would be too general, and the program would run in an almost endless loop
    jess = pyjess.Jess(templates) # Create a Jess instance and use it to query a molecule (a PDB structure) against the stored templates:
    query = jess.query(molecule=molecule, rmsd_threshold=rmsd_threshold, distance_cutoff=distance_cutoff, max_dynamic_distance=max_dynamic_distance, max_candidates=max_candidates, best_match=True) # query is pyjess.Query object which is an iterator over pyjess.Hits

    # A template is not encoded as coordinates, rather as a set of constraints. For example, it would not contain the exact positions of THR and ASN atoms, but instructions like "Cα of ASN should be X angstrom away from the Cα of THR plus the allowed distance."
    # Multiple solutions = Mathces to a template, satisfying all constraints may therefore exist
    # Jess produces matches to templates by looking for any combination of atoms, residue_types, elements etc. and ANY positions which satisfy the constraints in the template
    # thus the solutions that Jess finds are NOT identical to the template at all - rather they are all possible solutions to the set constraints. Solutions may completely differ from the template geometry or atom composition if allowed by the set constraints.
    # by setting best_match=True we turn on filtering by rmsd to return only the best match for every molecule template pair. Currently I dont expose best_match to the user and assume it as default. This should be the only use case (I think)
    # Thus we hope to return the one solution to the constraints which most closely resembles the original template - this is not guaranteed of course

    matches: List[Match] = []
    for hit in query: # hit is pyjess.Hit
        template = id_to_template[hit.template.id]
        matches.append(Match(hit=hit, template=template))

    def _check_completeness(matches: List[Match]):
        # only after all templates of a certain size have been scanned could we compute the complete tag
        # This requries cluster and mcsa.id to be set! Otherwise I assume there is no cluster and therefore the match is complete by default!
        groupable_matches = []
        lone_matches = []

        for match in matches:
            if match.template.mcsa_id is not None and match.template.cluster is not None:
                groupable_matches.append(match)
            else:
                lone_matches.append(match)

        # Group hit objects by the triple (hit.template.m-csa, hit.template.cluster.id, hit.template.dimension)
        def get_key(obj: Match) -> Tuple[int, int, int]:
            return obj.template.mcsa_id, obj.template.cluster.id, obj.template.dimension # type: ignore

        grouped_matches = [list(g) for _, g in itertools.groupby(sorted(groupable_matches, key=get_key), get_key)]

        for cluster_matches in grouped_matches:
            # For each query check if all Templates assigned to the same cluster targeted that structure
            #
            # TODO report statistics on this: This percentage of queries had a complete active site as reported by the complete tag
            # Filter this by template clusters with >1 member of course or report seperately by the number of clustermembers
            # or say like: This template cluster was always complete while this template cluster was only complete X times out of Y Queries matched to one member
            #
            # check if all the cluster members up to and including cluster_size are present in the group,
            indexed_possible_cluster_members = list(range(cluster_matches[0].template.cluster.size)) # type: ignore
            possible_cluster_members = [x+1 for x in indexed_possible_cluster_members]
            
            found_cluster_members = [match.template.cluster.member for match in cluster_matches] # type: ignore
            found_cluster_members.sort()

            if found_cluster_members == possible_cluster_members:
                for match in cluster_matches:
                    match.complete = True

        for match in lone_matches:
            match.complete=True
    
    _check_completeness(matches)

    return matches

class Matcher:
    def __init__(self, templates: List[AnnotatedTemplate], jess_params: Dict[int, Dict[str, float]] = {}, conservation_cutoff: int = 0, warn: bool = False, verbose: bool = False, skip_smaller_hits: bool = False, match_small_templates: bool = False, cpus: int = 0):

        self.templates = templates
        self.cpus = cpus
        self.conservation_cutoff = conservation_cutoff
        self.warn = warn
        self.verbose = verbose
        self.skip_smaller_hits = skip_smaller_hits
        self.match_small_templates = match_small_templates

        if self.cpus == 0:
            self.cpus = -2 + self.cpus + (os.cpu_count() or 1)
        if self.cpus < 0:
            self.cpus = (os.cpu_count() or 1)

        self.verbose_print(f'PyJess Version: {pyjess.__version__}')
        self.verbose_print(f'Running on {self.cpus} Thread(s)')
        self.verbose_print(f'Warnings are set to {self.warn}')
        self.verbose_print(f'Skip_smaller_hits search is set to {self.skip_smaller_hits}')

        if self.conservation_cutoff:
            self.verbose_print(f'Conservation Cutoff set to {self.conservation_cutoff}')

        if jess_params:
                self.jess_params = jess_params
        else:
            ####### Default Jess parameters by Size ##################################
            self.jess_params = {
                3: {'rmsd': 2, 'distance': 0.9, 'max_dynamic_distance': 0.9},
                4: {'rmsd': 2, 'distance': 1.7, 'max_dynamic_distance': 1.7},
                5: {'rmsd': 2, 'distance': 2.0, 'max_dynamic_distance': 2.0},
                6: {'rmsd': 2, 'distance': 2.0, 'max_dynamic_distance': 2.0},
                7: {'rmsd': 2, 'distance': 2.0, 'max_dynamic_distance': 2.0},
                8: {'rmsd': 2, 'distance': 2.0, 'max_dynamic_distance': 2.0}}

        # check each template and if it passes add it to the dictionary of templates
        self.templates_by_effective_size: Dict[int, List[AnnotatedTemplate]] = collections.defaultdict(list) # Dictionary of List of AnnoatedTemplate objects grouped by Template.effective_size as keys
        for template in templates:
            if check_template(template, warn=self.warn): # returns True if the Template passed all checks or if warn is set to False
                self.templates_by_effective_size[template.effective_size].append(template)

        if self.verbose:
            template_number_dict: Dict[int, int] = {}
            for size, template_list in self.templates_by_effective_size.items():
                if not self.match_small_templates and size < 3:
                    continue
                template_number_dict[size] = len(template_list)
            print(f'Templates by effective size: {collections.OrderedDict(sorted(template_number_dict.items()))}')
        
        self.template_effective_sizes = list(self.templates_by_effective_size.keys())
        self.template_effective_sizes.sort(reverse=True) # get a list of template_sizes in decending order

        if self.warn:
            smaller_sizes = [i for i in self.template_effective_sizes if i < 3]
            if smaller_sizes:
                small_templates = []
                for i in smaller_sizes:
                    small_templates.extend(self.templates_by_effective_size[i])

                if self.match_small_templates:
                    warnings.warn(f'{len(small_templates)} Templates with an effective size smaller than 3 defined sidechain residues were supplied.\nFor small templates Jess parameters for templates of 3 residues will be used.')
                else:
                    warnings.warn(f'{len(small_templates)} Templates with an effective size smaller than 3 defined sidechain residues were supplied.\nThese will be excluded since these templates are too general.')

                self.verbose_print('The templates with the following ids are too small:')
                self.verbose_print(small_templates)

    def verbose_print(self, *args):
        if self.verbose:
            print(*args)

    def run(self, molecules: List[pyjess.Molecule]) -> Dict[pyjess.Molecule, List[Match]]:
        processed_molecules = collections.defaultdict(list)

        with ThreadPool(self.cpus) as pool:

            for template_size in self.template_effective_sizes:
                templates = self.templates_by_effective_size[template_size]
                
                if template_size < 3:
                    if not self.match_small_templates:
                        continue # skip the rest of the for-loop
                    parameter_size = 3

                elif template_size > 8:
                    parameter_size = 8

                else:
                    parameter_size = template_size

                ################# Running Jess ###################################
                
                self.verbose_print(f'Now matching query structure(s) to template of size {template_size}')
                self.verbose_print(f'jess parameters are: {self.jess_params[parameter_size]}')

                rmsd = self.jess_params[parameter_size]['rmsd']
                distance = self.jess_params[parameter_size]['distance']
                max_dynamic_distance = self.jess_params[parameter_size]['max_dynamic_distance']
                
                _single_query_run = functools.partial(single_query_run, templates=templates, rmsd_threshold=rmsd, distance_cutoff=distance, max_dynamic_distance=max_dynamic_distance)

                if self.skip_smaller_hits:
                    query_molecules = [ molecule for molecule in molecules if molecule not in processed_molecules ]
                else:
                    query_molecules = molecules

                total_matches = 0 # counter for all matches found over all molecules passed 
                all_matches = pool.map(_single_query_run, query_molecules) # all_matches are a list of Match objects
                for molecule, matches in zip(query_molecules, all_matches):
                    processed_molecules[molecule].extend(matches)
                    total_matches += len(matches)

                self.verbose_print(f'{total_matches} matches were found for templates of effective size {template_size}')

        return processed_molecules

def main(argv: Optional[List[str]] = None, stderr=sys.stderr):

    class ReadListAction(argparse.Action):
        def __call__(self, parser, namespace, values, option_string=None):
            with values.open('r') as f:
                for line in f:
                    dest = getattr(namespace, self.dest)
                    dest.append(Path(line.strip()))

    parser = argparse.ArgumentParser(
        formatter_class=argparse.ArgumentDefaultsHelpFormatter,
        exit_on_error=False
    )
    # positional arguments
    # parser.add_argument('input', type=Path, help='File path of a single query pdb file OR File path of a file listing multiple query pdb files seperated by linebreaks')
    parser.add_argument('-o', "--output", required=True, type=Path, help='Output tsv file to which results should get written')
    parser.add_argument('--pdbs', type=Path, help='Output directory to which results should get written', default=None)

    # inputs: either a list of paths, or directly a path (or any combination)
    parser.add_argument('-i', '--input', type=Path, help='File path to a PDB file to use as query', action="append", dest="files", default=[])
    parser.add_argument('-l', '--list', type=Path, help="File containing a list of PDB files to read", action=ReadListAction, dest="files")

    # optional arguments
    parser.add_argument('-j', '--jess', nargs = 3, default=None, type=float, help='Fixed Jess parameters for all templates. Jess space seperated parameters rmsd, distance, max_dynamic_distance')
    parser.add_argument('-t', '--template-dir', type=Path, default=None, help='Path to directory containing jess templates. This directory will be recursively searched.')
    parser.add_argument('-v', '--verbose', default=False, action="store_true", help='If process information and time progress should be printed to the command line')
    parser.add_argument('-w', '--warn', default=False, action="store_true", help='If warings about bad template processing or suspicous and missing annotations should be raised')
    parser.add_argument('-q', '--include-query', default=False, action="store_true", help='Include the query structure together with the hits in the pdb output')
    parser.add_argument('-c', '--conservation-cutoff', type=float, default=0, help='Atoms with a value in the B-factor column below this cutoff will be excluded form matching to the templates. Useful for predicted structures.')
    parser.add_argument('-n', '--n-jobs', type=int, default=0, help='The number of threads to run in parallel. Pass 1 to run everything in the main thread, 0 to automatically select a suitable number, or any postive number. Negative numbers take all.')

    parser.add_argument('--transform', default=False, action="store_true", help='Transform the coordinate system of the hits to that of the template in the pdb output')
    parser.add_argument('--skip-smaller-hits', default=False, action="store_true", help='If True, do not search with smaller templates if larger templates have already found hits.')
    parser.add_argument('--match-small-templates', default=False, action="store_true", help='If true, templates with less then 3 defined sidechain residues will still be matched.')
    
    args = parser.parse_args(args=argv)
    
    if not args.files:
        parser.error("No input files were passed. Use the -i and/or -l flags to pass input.")

    jess_params = {}
    if args.jess:
        jess_params_list = [i for i in args.jess]
        # jess parameters
        # we use different parameters for different template residue numbers - higher number more generous parameters
        rmsd = jess_params_list[0] # in Angstrom, typcically set to 2
        distance = jess_params_list[1] # in Angstrom between 1.0 and 1.5 - lower is more strict. This changes with template size
        max_dynamic_distance = jess_params_list[2] # if equal to distance dynamic is off: this option is currenlty dysfunctional
    
        jess_params = {
            3: {'rmsd': rmsd, 'distance': distance, 'max_dynamic_distance': max_dynamic_distance},
            4: {'rmsd': rmsd, 'distance': distance, 'max_dynamic_distance': max_dynamic_distance},
            5: {'rmsd': rmsd, 'distance': distance, 'max_dynamic_distance': max_dynamic_distance},
            6: {'rmsd': rmsd, 'distance': distance, 'max_dynamic_distance': max_dynamic_distance},
            7: {'rmsd': rmsd, 'distance': distance, 'max_dynamic_distance': max_dynamic_distance},
            8: {'rmsd': rmsd, 'distance': distance, 'max_dynamic_distance': max_dynamic_distance}}

    try:

        ######## Loading all query molecules #############################
        molecules: List[pyjess.Molecule] = []
        stem_counter: Dict[str, int] = collections.defaultdict(int)
        for molecule_path in args.files:
            stem = Path(molecule_path).stem
            stem_counter[stem] += 1
            if stem_counter[stem] > 1:
                # In case the same stem occurs multiple times, create a unique ID using the stem and a running number starting from 2
                unique_id = f"{stem}_{stem_counter[stem]}"
            else:
                unique_id = stem

            molecule = pyjess.Molecule.load(str(molecule_path), id=unique_id) # by default it will stop at ENDMDL
            if args.conservation_cutoff: # if the pyjess cutoff function changes to stop it from making a copy if conservation cutoff is 0 or false or something I can drop this
                molecule = molecule.conserved(args.conservation_cutoff)
                # conserved is a method called on a molecule object that returns a filtered molecule
                # atoms with a temperature-factor BELOW the conservation cutoff will be excluded
            if molecule:
                molecules.append(molecule) # load a molecule and filter it by conservation_cutoff
            elif args.warn:
                warnings.warn(f"received an empty molecule from {molecule_path}")

        if not molecules and args.warn:
            warnings.warn("received no molecules from input")

        if args.template_dir:
            templates = list(load_templates(template_dir=Path(args.template_dir), warn=args.warn, verbose=args.verbose))
        else:
            templates = list(load_templates(warn=args.warn, verbose=args.verbose))

        ############ Initialize Matcher object ################################
        matcher = Matcher(templates=templates, jess_params=jess_params, conservation_cutoff=args.conservation_cutoff, warn=args.warn, verbose=args.verbose, skip_smaller_hits=args.skip_smaller_hits, match_small_templates=args.match_small_templates, cpus=args.n_jobs)

        ############ Call Matcher.run ##########################################
        processed_molecules = matcher.run(molecules = molecules)

        ######### Writing Output ##########################################
        out_tsv = args.output

        if not out_tsv.parent.exists():
            out_tsv.parent.mkdir(parents=True, exist_ok=True)
            if args.warn:
                warnings.warn(f"{out_tsv.parent.resolve()} directory tree to output tsv file did not exist and was created")
        elif out_tsv.exists() and args.warn:
            warnings.warn(f'The specified output tsv file {out_tsv.resolve()} already exists and will be overwritten!')

        if args.verbose:
            print(f'Writing output to {out_tsv.resolve()}')

        with open(out_tsv, 'w', newline='', encoding ="utf-8") as tsvfile:
            for index, (molecule, matches) in enumerate(processed_molecules.items()):
                for jndex, match in enumerate(matches):
                    i = index+jndex
                    match.index = jndex+1 # 1 indexed matches per query
                    match.dump(tsvfile, header=i==0) # one line per match, write header only for the first match

        def write_hits2_pdb(matches: List[Match], filename: str, outdir: Path):
            outdir.mkdir(parents=True, exist_ok=True)
            # make sure molecule.id is unique!
            with open(Path(outdir, f'{filename}_matches.pdb'), 'w', encoding ="utf-8") as pdbfile:
                if args.include_query: # write the molecule structure to the top of the pdb output too
                    for i, match in enumerate(matches):
                        match.dump2pdb(pdbfile, include_query=i==0, transform=args.transform)
                else:
                    for i, match in enumerate(matches):
                        match.dump2pdb(pdbfile, transform=args.transform)

        if args.pdbs:
            for molecule, matches in processed_molecules.items():
                write_hits2_pdb(matches=matches, filename=molecule.id, outdir=args.pdbs) # type: ignore

    except IsADirectoryError as exc:
        print("File is a directory:", exc.filename, file=stderr)
        return exc.errno

    except FileNotFoundError as exc:
        print("Failed to find file:", exc.filename, file=stderr)
        return exc.errno

    except FileExistsError as exc:
        print("File already exists:", exc.filename, file=stderr)
        return exc.errno

    return 0

if __name__ == "__main__":
    sys.exit(main())