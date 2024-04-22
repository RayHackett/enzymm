import argparse
from pathlib import Path
from importlib.resources import files
from typing import List, Set, Tuple, Dict, Optional, Union, TextIO, Iterator, Iterable, IO
from functools import cached_property
from dataclasses import dataclass, field
import tempfile
import warnings
import csv
import numpy as np
import itertools
import io

import pyjess # type: ignore # cython port of Jess to python package

from .jess.filter_hits import filter_hits
from .template import Template, Vec3, load_templates, check_template, group_templates_by_size
from .utils import chunks, ranked_argsort

__all__ = [
    "Match"
]

@dataclass
class Match:
    """"Class for storing annotated Jess hits. Wrapper around pyjess.Hit and the original Template object"""
    hit: pyjess.Hit
    template: Template
    complete: bool = field(default=False) # initialize as false but change it later
    # maybe have the raw PDB output for use with pymol too

    def dumps(self) -> str:
        buffer = io.StringIO()
        self.dump(buffer)
        return buffer.getvalue() # returns entire content temporary file object as a string

    def dump2pdb(self, file: IO[str]):
        file.write(f'REMARK TEMPLATE_ID {self.template.id}\n')
        file.write(f'REMARK MOLECULE_ID {self.hit.molecule.id}\n')

        one_char_elements = {'H', 'B', 'C', 'N', 'O', 'F', 'P', 'S', 'K', 'V', 'Y', 'I', 'W', 'U'}
        for atom in self.hit.atoms(transform=False): # keep the atoms in the query reference frame!
            if atom.element in one_char_elements:
                file.write(f"ATOM  {atom.serial:>5}  {atom.name:<3s}{atom.altloc if atom.altloc is not None else '':<1}{atom.residue_name:<3}{atom.chain_id:>2}{atom.residue_number:>4}{atom.insertion_code:1s}   {atom.x:>8.3f}{atom.y:>8.3f}{atom.z:>8.3f}{atom.occupancy:>6.2f}{atom.temperature_factor:>6.2f}      {atom.segment:<4s}{atom.element:>2s} \n")
            else:
                file.write(f"ATOM  {atom.serial:>5} {atom.name:<4s}{atom.altloc if atom.altloc is not None else '':<1}{atom.residue_name:<3}{atom.chain_id:>2}{atom.residue_number:>4}{atom.insertion_code:1s}   {atom.x:>8.3f}{atom.y:>8.3f}{atom.z:>8.3f}{atom.occupancy:>6.2f}{atom.temperature_factor:>6.2f}      {atom.segment:<4s}{atom.element:>2s} \n")
        file.write('END\n\n')

    def dump(self, file: IO[str], header: bool = False):
        writer = csv.writer(file, dialect="excel-tab", delimiter='\t', lineterminator='\n')
        if header:
                writer.writerow([
                    "query_id",
                    "template_id",
                    "template_size",
                    "template_true_size",
                    "template_mcsa_id",
                    "template_uniprot_id",
                    "template_pdb_id",
                    "template_ec",
                    "template_cath",
                    "template_multimeric",
                    "query_multimeric",
                    "rmsd",
                    "log_evalue",
                    "orientation",
                    "preserved_order",
                    "completeness",
                    "matched_residues"])

        writer.writerow([
                str(self.hit.molecule.id),
                str(self.template.id),
                str(self.template.size),
                str(self.template.true_size),
                str(self.template.mcsa_id),
                str(self.template.uniprot_id),
                str(self.template.pdb_id),
                ",".join(self.template.ec if self.template.ec is not None else ''),
                ",".join(self.template.cath if self.template.cath else ''),
                str(self.template.multimeric),
                str(self.multimeric),
                str(self.hit.rmsd),
                str(self.hit.log_evalue),
                str(self.orientation),
                str(self.preserved_resid_order),
                str(self.complete),
                (','.join('_'.join(t) for t in self.matched_residues))])

    @cached_property
    def atom_triplets(self) -> List[List[pyjess.Atom]]:
        # list with matched residues
        # # Hit.atoms is a list of matched atoms with all info on residue numbers and residue chain ids and atom types, this should conserve order if Hit.atoms is a list!!!
        atom_triplets : List[List[pyjess.Atom]] = []
        for atom_triplet in chunks(self.hit.atoms(), 3): # yield chunks of 3 atoms each
            # check if all three atoms belong to the same residue by adding a tuple of their residue defining properties to a set
            unique_residues = { (atom.residue_name, atom.chain_id, atom.residue_number) for atom in atom_triplet }
            if len(unique_residues) != 1:
                raise ValueError('Mixed up atom triplets. The atoms come from different residues!')
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
        if self.multimeric:
            return False
        else:
            # Now extract relative atom order in hit
            return ranked_argsort([atom_triplet[0].residue_number for atom_triplet in self.atom_triplets]) == self.template.relative_order

    @cached_property
    def match_vector_list(cls) -> List[Vec3]:
        vector_list = []
        for residue_index, residue in enumerate(cls.template.residues):
            first_atom_index, second_atom_index = residue.orientation_vector_indices
            if second_atom_index == 9:
                middle_atom = cls.atom_triplets[residue_index][first_atom_index]
                side1, side2 = [atom for atom in cls.atom_triplets[residue_index] if atom != middle_atom]
                midpoint = [(side1.x + side2.x) / 2, (side1.y + side2.y) / 2, (side1.z + side2.z) / 2]
                vector_list.append(Vec3(midpoint[0]- middle_atom.x, midpoint[1] - middle_atom.y, midpoint[2]- middle_atom.z))
            else:
                # Calculate orientation vector going from first_atom to second_atom_index
                first_atom = cls.atom_triplets[residue_index][first_atom_index]
                second_atom = cls.atom_triplets[residue_index][second_atom_index]
                vector_list.append(Vec3(second_atom.x - first_atom.x, second_atom.y - second_atom.y, second_atom.z - first_atom.z))
        return vector_list

    @cached_property
    def template_vector_list(self) -> List[Vec3]:
        return [res.orientation_vector for res in self.template.residues]
    
    @property
    def orientation(self) -> float: # average angle
        if len(self.template_vector_list) != len(self.match_vector_list):
            raise ValueError('Vector lists for Template and matching Query structure had different lengths.')

        def angle_between(v1: Vec3, v2: Vec3): # from https://stackoverflow.com/questions/2827393/angles-between-two-n-dimensional-vectors-in-python
            """ Returns the angle in radians between vectors 'v1' and 'v2':
            """
            try:
                a1 = np.array([v1.x, v1.y, v1.z])
                a2 = np.array([v2.x, v2.y, v2.z])
                a1_u = a1 / np.linalg.norm(a1) # unit vector
                a2_u = a2 / np.linalg.norm(a2) # unit vector
                return np.arccos(np.clip(a1_u @ a2_u, -1.0, 1.0))
            except Exception as exc:
                raise ValueError(f'vectorx {v1.x}, vectory {v1.y}, vectorz {v1.z}') from exc

        # now calculate the angle between the vector of the template and the query per residue
        angle_list = []
        for i in range(len(self.template_vector_list)):
            angle_list.append(angle_between(self.match_vector_list[i], self.template_vector_list[i]))

        return np.mean(angle_list)

def _single_query_run(molecule: pyjess.Molecule, pyjess_templates: Iterable[pyjess.Template], id_to_template: Dict[str, Template], rmsd: float = 2.0, distance: float = 3.0, max_dynamic_distance: float = 3.0, max_candidates: int = 10000):
    # killswitch is controlled by max_candidates. Internal default is currently 1000
    # killswitch serves to limit the iterations in cases where the template would be too general, and the program would run in an almost endless loop
    jess = pyjess.Jess(pyjess_templates) # Create a Jess instance and use it to query a molecule (a PDB structure) against the stored templates:
    query = jess.query(molecule, rmsd_threshold=rmsd, distance_cutoff=distance, max_dynamic_distance=max_dynamic_distance, max_candidates=max_candidates)

    hits = list(query)
    # A template is not encoded as coordinates, rather as a set of constraints. For example, it would not contain the exact positions of THR and ASN atoms, but instructions like "Cα of ASN should be X angstrom away from the Cα of THR plus the allowed distance."
    # Multiple solutions = Mathces to a template, satisfying all constraints may therefore exist
    # Jess produces matches to templates by looking for any combination of atoms, residue_types, elements etc. and ANY positions which satisfy the constraints in the template
    # thus the solutions that Jess finds are NOT identical to the template at all - rather they are all possible solutions to the set constraints. Solutions may completely differ from the template geometry or atom composition if allowed by the set constraints.
    # by filtering for e-value we hope to return the solution to the constraints which most closely resembles the original template.
    matches: List[Match] = []
    if hits: # if any hits were found
        # best hits is a list of pyjess.Hit objects
        best_hits =  filter_hits(hits) # retain only the hit with the lowest e-value for each query-template pair

        if best_hits:
            for hit in best_hits:
                template = id_to_template[hit.template.id]
                matches.append(Match(hit=hit, template=template))

        def _check_completeness(matches: List[Match]):
            # only after all templates of a certain size have been scanned could we compute the complete tag
            
            # Group hit objects by the pair (hit.template.m-csa, first digit of hit.template.cluster)
            def get_key(obj: Match) -> Tuple[int, int]:
                return obj.template.mcsa_id, obj.template.cluster_id # split on _ and take the first digit which is the cluster id

            grouped_matches = [list(g) for _, g in itertools.groupby(sorted(matches, key=get_key), get_key)]

            for cluster_matches in grouped_matches:
                # For each query check if all Templates assigned to the same cluster targeted that structure
                #
                # TODO report statistics on this: This percentage of queries had a complete active site as reported by the complete tag
                # Filter this by template clusters with >1 member of course or report seperately by the number of clustermembers
                # or say like: This template cluster was always complete while this template cluster was only complete X times out of Y Queries matched to one member
                #
                # check if all the cluster members up to and including cluster_size are present in the group,
                indexed_possible_cluster_members = list(range(cluster_matches[0].template.cluster_size))
                possible_cluster_members = [x+1 for x in indexed_possible_cluster_members]
                
                found_cluster_members = [match.template.cluster_member for match in cluster_matches] # second number is the current cluster member number

                if found_cluster_members == possible_cluster_members:
                    for match in cluster_matches:
                        match.complete = True
        
        _check_completeness(matches)

    return matches

def write_matches_to_tsv(matches: List[Match], filename: str, outdir: Path):
    results_path = Path(outdir, "results/")
    results_path.mkdir(parents=True, exist_ok=True)
    # TODO somehow we need to give unique filenames in case molecule.id is not unique
    with open(Path(results_path, f'{filename}_matches.tsv'), 'w', newline='', encoding ="utf-8") as tsvfile:
        for i, match in enumerate(matches):
            match.dump(tsvfile, header=i==0) # one line per match, write header only for the first match too

def write_hits2_pdb(matches: List[Match], filename: str, outdir: Path):
    results_path = Path(outdir, "results/")
    results_path.mkdir(parents=True, exist_ok=True)
    # TODO somehow we need to give unique filenames in case molecule.id is not unique
    with open(Path(results_path, f'{filename}_matches.pdb'), 'w', encoding ="utf-8") as pdbfile:
        for i, match in enumerate(matches):
            match.dump2pdb(pdbfile)

def matcher_run(query_path: Path, template_path: Path, jess_params: Dict[int, Dict[str, float]], outdir: Path, conservation_cutoff: int = 0, warn: bool = True, verbose: bool = False, skip_smaller_hits: bool = False):

    if verbose:
        print(f'Warnings are set to {warn}')
        print(f'Skip_smaller_hits search is set to {skip_smaller_hits}')

    ####### Checking outdir ##########################################
    try:
        if not outdir.is_dir():
            outdir.mkdir(parents=True, exist_ok=True)
            if warn:
                warnings.warn(f"'{outdir.resolve()}' directory did not exist and was created")
        if verbose:
            print(f'Writing output to {Path(outdir).resolve()}')

    except ValueError as exc: # permission error or something that did not look like a path was passed?
        raise ValueError(f'{outdir} passed to output was not a valid directory path') from exc
    
    ######## Checking query input files #################################
    def pdb_file_check(path_to_check: Path) -> bool:
        return Path(path_to_check).exists() and Path(path_to_check).suffix.lower() in {".pdb", ".ent"}

    molecule_paths: Set[Path] = set()
    if pdb_file_check(query_path):
        molecule_paths.add(query_path)
    else:
        try:
            with open(Path(query_path), 'r') as f:
                lines = f.readlines()
            for line in lines:
                if pdb_file_check(Path(line.strip())):
                    molecule_paths.add(Path(line.strip()))
        except:
            raise FileNotFoundError(f'File {query_path} did not exist or did not contain any paths to files with the expected .pdb or .ent extensions')

    ######## Loading all query molecules #############################
    molecules: Set[pyjess.Molecule] = set()
    if conservation_cutoff:
        if verbose:
            print(f'Conservation Cutoff set to {conservation_cutoff}')
        for molecule_path in molecule_paths:
            # conserved is a method called on a molecule object that returns a filtered molecule
            # atoms with a temperature-factor BELOW the conservation cutoff will be excluded
            molecules.add(pyjess.Molecule.load(molecule_path).conserved(conservation_cutoff)) # load a molecule and filter it by conservation_cutoff
    else:
        for molecule_path in molecule_paths:
            molecules.add(pyjess.Molecule.load(molecule_path)) # load a molecule

    ######## Checking template path ##################################

    if template_path:
        if verbose:
            print(f'loading supplied template files from {template_path.resolve()}')
        template_tuples = list(load_templates(template_path, warn=warn))
    else:
        if verbose:
            print('loading default template files')
        template_tuples = list(load_templates(warn=warn)) # default templates

    # check each template and if it passes add it to the mapping dictionary and transform it to a pyjess.Template and add it to a list of pyjess.Template objects
    template_size_to_pyjess_template_id: Dict[int, List[pyjess.Template]] = {} # Dictionary of List of pyjess.Template objects grouped by Template.size as keys
    id_to_template: Dict[str, Template] = {} # mapping pyjess.Template ID to my Template objects
    for i, template_tuple in enumerate(template_tuples):
        if check_template(template_tuple, warn=warn): # returns True if the Template passed all checks or if warn is set to False
            pyjess_template = template_tuple[0].to_pyjess_template()
            id_to_template[pyjess_template.id] = template_tuple[0]
            if template_tuple[0].size not in template_size_to_pyjess_template_id:
                template_size_to_pyjess_template_id[template_tuple[0].size] = []
            template_size_to_pyjess_template_id[template_tuple[0].size].append(pyjess_template)
    
    template_sizes = list(template_size_to_pyjess_template_id.keys())
    template_sizes.sort(reverse=True) # get a list of template_sizes in decending order

    if warn:
        for i in template_sizes:
            if i < 3:
                print('Templates with a size smaller than 3 defined, sidechain residues were supplied. These will be excluded since these templates are too general')

    remaining_molecules = molecules
    processed_molecules: Dict[pyjess.Molecule, List[Match]] = dict()
    for template_size in template_sizes:
        pyjess_templates = template_size_to_pyjess_template_id[template_size]
        
        if template_size < 3:
            if warn:
                print('The following templates are too small and will be skipped:')
                for tmp in pyjess_templates:
                    print(tmp.id)
            continue # skip the rest of the for-loop

        elif template_size > 8:
            if warn:
                print('Templates with more than 8 residues were passed. Jess parameters for templates of 8 residues will be used.')
            parameter_size = 8

        else:
            parameter_size = template_size

        ################# Running Jess ###################################
        
        if verbose:
            print(f'Now matching query structure(s) to template of size {template_size}')
            print(f'jess parameters are: {jess_params[parameter_size]}')

        rmsd = jess_params[parameter_size]['rmsd']
        distance = jess_params[parameter_size]['distance']
        max_dynamic_distance = jess_params[parameter_size]['max_dynamic_distance']

        total_matches = 0 # counter for all matches found over all molecules passed
        # iterate over all the query molecules passed
        for molecule in remaining_molecules:
            matches = _single_query_run(molecule=molecule, pyjess_templates=pyjess_templates, id_to_template=id_to_template, rmsd=rmsd, distance=distance, max_dynamic_distance=max_dynamic_distance)
            if matches:
                if molecule not in processed_molecules:
                    processed_molecules[molecule] = []
                processed_molecules[molecule].extend(matches)
                total_matches += len(matches)

        if skip_smaller_hits: # if skip_smaller_hits is true, then do not search with smaller templates if hits with larger ones have been found.
            remaining_molecules = remaining_molecules.difference(set(processed_molecules.keys()))

        if verbose:
            print(f"found {total_matches} matches with templates of size {template_size}")

    for molecule, matches in processed_molecules.items():
        write_matches_to_tsv(matches=matches, filename=molecule.id, outdir=outdir)
        write_hits2_pdb(matches=matches, filename=molecule.id, outdir=outdir)

    # TODO what should i retrun? a path to the output file? a bool if any hits were found?
    # return
        
if __name__ == "__main__":
    
    parser = argparse.ArgumentParser()
    # positional arguments
    parser.add_argument('input', type=Path, help='File path of a single query pdb file OR File path of a file listing multiple query pdb files seperated by linebreaks')
    parser.add_argument('output_dir', type=Path, help='Output Directory to which results should get written')

    # optional arguments
    parser.add_argument('-j', '--jess', nargs = 4, default=None, type=float, help='Fixed Jess parameters for all templates. Jess space seperated parameters rmsd, distance, max_dynamic_distance, score_cutoff')
    parser.add_argument('-t', '--template_dir', type=Path, default=None, help='Path to directory containing jess templates. This directory will be recursively searched.')
    parser.add_argument('-s', '--skip_smaller_hits', default=False, action="store_true", help='If True, do not search with smaller templates if larger templates have already found hits.')
    parser.add_argument('-v', '--verbose', default=False, action="store_true", help='If process information and time progress should be printed to the command line')
    parser.add_argument('-w', '--warn', default=False, action="store_true", help='If warings about bad template processing or suspicous and missing annotations should be raised')
    parser.add_argument('--conservation_cutoff', default=None, help='Atoms with a value in the B-factor column below this cutoff will be excluded form matching to the templates')
    args = parser.parse_args()
    
    query_path = args.input
    template_path = args.template_dir
    conservation_cutoff =  args.conservation_cutoff # Reads B-factor column in .pdb files: atoms below this cutoff will be disregarded. Could be pLDDT for example

    if args.jess:
        jess_params_list = [i for i in args.jess]
        # jess parameters
        # we use different parameters for different template residue numbers - higher number more generous parameters
        rmsd = jess_params_list[0] # in Angstrom, typcically set to 2
        distance = jess_params_list[1] # in Angstrom between 1.0 and 1.5 - lower is more strict. This changes with template size
        max_dynamic_distance = jess_params_list[2] # if equal to distance dynamic is off: this option is currenlty dysfunctional
        

        jess_params = {
            3: {'rmsd': rmsd, 'distance': distance, 'max_dynamic_distance': distance},
            4: {'rmsd': rmsd, 'distance': distance, 'max_dynamic_distance': distance},
            5: {'rmsd': rmsd, 'distance': distance, 'max_dynamic_distance': distance},
            6: {'rmsd': rmsd, 'distance': distance, 'max_dynamic_distance': distance},
            7: {'rmsd': rmsd, 'distance': distance, 'max_dynamic_distance': distance},
            8: {'rmsd': rmsd, 'distance': distance, 'max_dynamic_distance': distance}}

    else:
        ####### Default Jess parameters by Size ##################################
        jess_params = {
            3: {'rmsd': 2, 'distance': 1, 'max_dynamic_distance': 1},
            4: {'rmsd': 2, 'distance': 1.5, 'max_dynamic_distance': 1.5},
            5: {'rmsd': 2, 'distance': 1.5, 'max_dynamic_distance': 1.5},
            6: {'rmsd': 2, 'distance': 1.5, 'max_dynamic_distance': 1.5},
            7: {'rmsd': 2, 'distance': 1.5, 'max_dynamic_distance': 1.5},
            8: {'rmsd': 2, 'distance': 1.5, 'max_dynamic_distance': 1.5}}

    skip_smaller_hits = args.skip_smaller_hits
    outdir = args.output_dir
    warn = args.warn
    verbose = args.verbose

    matcher_run(query_path=query_path, template_path=template_path, jess_params=jess_params, outdir=outdir, conservation_cutoff=conservation_cutoff, warn=warn, verbose=verbose, skip_smaller_hits=skip_smaller_hits)