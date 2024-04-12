import argparse
from pathlib import Path
from importlib.resources import files
from typing import List, Set, Tuple, Dict, Optional, Union, TextIO, Iterator, Iterable
from functools import cached_property
from dataclasses import dataclass, field
import tempfile
import warnings
import csv
import numpy as np

import pyjess # cython port of Jess to python package

from .jess.filter_hits import filter_hits
from .template import Template, Vec3, load_templates, check_template
from .utils import chunks, ranked_argsort

__all__ = [
    "Match"
]

@dataclass
class Match:
    """"Class for storing annotated Jess hits. Wrapper around pyjess.Hit and the original Template object"""
    hit: pyjess.Hit
    template: Template
    completeness: bool = field(default=False) # initialize as false but change it later
    # maybe have the raw PDB output for use with pymol too

    def dump(self) -> List[str]:
        return [str(self.hit.molecule.id),
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
                str(self.completeness),
                (','.join(','.join(t) for t in self.matched_residues))]

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

    @cached_property
    def matched_residues(self) -> Set[Tuple[str, str, str]]:
        return { (atom_triplet[0].residue_name, atom_triplet[0].chain_id, str(atom_triplet[0].residue_number)) for atom_triplet in self.atom_triplets }

    @cached_property
    def multimeric(self) -> bool:
        """`bool`: Boolean if the atoms in the hit stem from multiple protein chains
        """
        # note that these are pyjess atom objects!
        return all(atom.chain_id == self.hit.atoms()[0].chain_id for atom in self.hit.atoms())
        # TODO pyjess atoms need to be iterable

    @cached_property
    def preserved_resid_order(self) -> bool:
        """`bool`: Boolean if the residues in the template and in the matched query structure have the same relative order.
        This is a good filtering parameter but excludes hits on examples of convergent evolution or circular permutations
        """
        if self.multimeric:
            return False
        else:
            # Now extract relative atom order in hit
            return ranked_argsort([atom.residue_number for atom in self.hit.atoms()]) == Template.relative_order

    @cached_property
    def match_vector_list(cls) -> List[Vec3]:
        # dictionary in which the vectors from start to finish are defined for each aminoacid type
        # the orientation vector is calculated differently for different aminoacid types
        vector_atom_type_dict = {'GLY': ('C','O'),
                            'ALA': ('CA','CB'),
                            'VAL': ('CA','CB'),
                            'LEU': ('CA','CB'),
                            'ILE': ('CA','CB'),
                            'MET': ('CG','SD'),
                            'PHE': ('CZ','mid'),
                            'TYR': ('CZ','OH'),
                            'TRP': ('CZ2','NE1'),
                            'CYS': ('CB','SG'),
                            'PRO': ('C','O'),
                            'SER': ('CB','OG'),
                            'THR': ('CB','OG1'),
                            'ASN': ('CG','OD1'),
                            'GLN': ('CD','OE1'),
                            'LYS': ('CE','NZ'),
                            'ARG': ('CZ','mid'),
                            'HIS': ('CG','ND1'),
                            'ASP': ('CG','mid'),
                            'GLU': ('CD','mid'),
                            'PTM': ('CA','CB'),
                            'ANY': ('C','O')}

        vector_list = []
        for atom_triplet in cls.atom_triplets:
            vectup = vector_atom_type_dict[atom_triplet[0].residue_name]
            # In residues with two identical atoms, the vector is calculated between the middle atom and the mid point between the identical pair
            if vectup[1] == 'mid':
                try:
                    middle = next(atom for atom in atom_triplet if atom.name == vectup[0])
                    side1, side2 = [atom for atom in atom_triplet if atom != middle]
                    midpoint = [(side1.x + side2.x) / 2, (side1.y + side2.y) / 2, (side1.z + side2.z) / 2]
                    vector_list.append(Vec3(middle.x - midpoint[0], middle.y - midpoint[1], middle.z - midpoint[2]))
                except StopIteration:
                    raise ValueError(f"Failed to find middle atom for amino-acid {atom_triplet[0].residue_name!r}") from None
                    
            else: # from first atom to second atom
                try:
                    first_atom = next(atom for atom in atom_triplet if atom.name == vectup[0])
                except StopIteration:
                    raise ValueError(f'Failed to find first atom for amino-acid {atom_triplet[0].residue_name!r}') from None
                try:
                    second_atom = next(atom for atom in atom_triplet if atom.name == vectup[1])
                except StopIteration:
                    raise ValueError(f'Failed to find second atom for amino-acid {atom_triplet[0].residue_name!r}') from None
                    
                vector_list.append(Vec3(first_atom.x - second_atom.x, first_atom.y - second_atom.y, first_atom.z - second_atom.z))
        return vector_list

    @cached_property
    def template_vector_list(self) -> List[Vec3]:
        return [res.orientation_vector for res in self.template.residues]
    
    @cached_property
    def orientation(self) -> float: # average angle
        if len(self.template_vector_list) != len(self.match_vector_list):
            raise ValueError('Vector lists for Template and matching Query structure had different lengths.')

        def angle_between(v1: Vec3, v2: Vec3): # from https://stackoverflow.com/questions/2827393/angles-between-two-n-dimensional-vectors-in-python
            """ Returns the angle in radians between vectors 'v1' and 'v2':
            """
            v1_u = v1 / np.linalg.norm(v1.x, v1.y, v1.z) # unit vector
            v2_u = v2 / np.linalg.norm(v2.x, v2.y, v2.z) # unit vector
            return np.arccos(np.clip(np.dot(v1_u, v2_u), -1.0, 1.0))

        # now calculate the angle between the vector of the template and the query per residue
        angle_list = []
        for i in range(len(self.template_vector_list)):
            angle_list.append(angle_between(self.match_vector_list[i], self.template_vector_list[i]))

        return np.mean(angle_list)

def _single_query_run(outdir: Path, molecule_path: Path, pyjess_templates: Iterable[pyjess.Template], id_to_template: Dict[str, Template], rmsd: float = 2.0, distance: float = 3.0, max_dynamic_distance: float = 3.0, conservation_cutoff: float = 100, max_candidates: int = 10000):
    # TODO should this function return anything?
    molecule = pyjess.Molecule.loads(molecule_path.open().read()) # This controls the B-factor column and removes anything smaller
    # killswitch is controlled by max_candidates. Internal default is currently 10000
    jess = pyjess.Jess(pyjess_templates) # Create a Jess instance and use it to query a molecule (a PDB structure) against the stored templates:
    query = jess.query(molecule, rmsd_threshold=rmsd, distance_cutoff=distance, max_dynamic_distance=max_dynamic_distance, max_candidates=max_candidates)

    if query: # if any hits were found
        # best hits is a list of pyjess.Hit objects
        best_hits =  filter_hits(list(query)) # retain only the hit with the lowest e-value for each query-template pair

    matches: List[Match] = []
    for hit in best_hits:
        template = id_to_template[hit.template.id]
        matches.append(Match(hit=hit, template=template))

    # TODO what kind of method is this? how to I add a required attribute to a class object afterwards?
    def _check_completeness(matches):
        # only after all templates of a certain size have been scanned could we compute the completeness tag
        
        # Group hit objects by the pair (hit.template.m-csa, first digit of hit.template.cluster)
        def get_key(obj) -> Tuple[str, str]:
            return obj.template.msca_id, obj.template.cluster.split('.')[0] # split on . and take the first digit which is the cluster number

        grouped_hits = [list(g) for _, g in itertools.groupby(sorted(matches, key=get_key), get_key)]

        for cluster_hits in grouped_hits:
            # For each query check if all Templates assigned to the same cluster targeted that structure
            #
            # TODO report statistics on this: This percentage of queries had a complete active site as reported by the completeness tag
            # Filter this by template clusters with >1 member of course or report seperately by the number of clustermembers
            # or say like: This template cluster was always complete while this template cluster was only complete X times out of Y Queries matched to one member
            #
            # get the total number templates from the first hit in each cluster which is the third cluster digit
            total_clustermembers = int(cluster_hits[0].template.cluster.split('.')[2])
            # check if all the cluster members with all 2nd digit identiers up to and including total_clustermembers are present in the group,
            indexed_possible_cluster_members = list(range(total_clustermembers))
            possible_cluster_members = [x+1 for x in indexed_possible_cluster_members]
            
            found_cluster_members = [int(hit.template.cluster.split('.')[1]) for hit in cluster_hits] # second number is the current cluster member number

            if found_cluster_members == possible_cluster_members:
                for hit in cluster_hits:
                    hit.complete = True
    
    _check_completeness(matches)

    results_path = Path(outdir, "results/")
    results_path.mkdir(parents=True, exist_ok=True)
    # TODO somehow we need to give unique filenames in case molecule_path.stem is not unique
    with open(Path(results_path, f'{molecule_path.stem}_matches.tsv'), 'w', newline='', encoding ="utf-8") as tsvfile:
        writer = csv.writer(tsvfile, delimiter='\t', lineterminator='\n')
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
                        "query_uniprot_id",
                        "query_pdb_id",
                        "query_ec",
                        "query_cath",
                        "query_multimeric",
                        "rmsd",
                        "log_evalue",
                        "orientation",
                        "preserved_order",
                        "completeness",
                        "matched_residues"])
        for match in matches:
            writer.writerow(match.dump()) # one line per match

def matcher_run(query_path: Path, rmsd: float, distance: float, max_dynamic_distance: float, conservation_cutoff: float, outdir: Path, warn: bool = True, verbose: bool = False):

    ####### Checking outdir ##########################################
    try:
        if not outdir.is_dir():
            outdir.mkdir(parents=True, exist_ok=True)
            if warn:
                warnings.warn(f"'{outdir}' directory did not exist and was created")
        if verbose:
            print(f'Writing output to {Path(outdir).resolve()}')

    except ValueError as exc: # permission error or something that did not look like a path was passed?
        raise ValueError(f'{outdir} passed to output was not a valid directory path') from exc
    
    ######## Checking query input files #################################
    def pdb_file_check(path_to_check: Path) -> bool:
        return Path(path_to_check).exists() and Path(path_to_check).suffix.lower() in {".pdb", ".ent"}

    if pdb_file_check(query_path):
        molecule_paths = [query_path]
    else:
        try:
            molecule_paths = []
            with open(Path(query_path), 'r') as f:
                lines = f.readlines()
            for line in lines:
                if pdb_file_check(Path(line.strip())):
                    molecule_paths.append(Path(line.strip()))
        except:
            raise FileNotFoundError(f'File {query_path} did not exist or did not contain any paths to files with the expected .pdb or .ent extensions')

    ################# Running Jess ###################################
    # TODO set level of verbosity and then print if conditions apply
    if verbose:
        print('jess parameters are: ', rmsd, distance, max_dynamic_distance, conservation_cutoff)
    
    # we iterate over templates starting with the ones with the largest ones
    for template_res_num in [8, 7, 6, 5, 4, 3]:
        if verbose:
            print('Now on template res num {template_res_num}')
        # load the templates with a given residue number, yields tuple with (Template,Path)
        template_tuples = list(load_templates(Path(str(files(__package__).joinpath('jess_templates_20230210', f'{template_res_num}_residues'))), warn=False))
        
        # check each template and if it passes add it to the mapping dictionary and transform it to a pyjess.Template and add it to a list of pyjess.Template objects
        pyjess_templates: List[pyjess.Template] = [] # List of pyjess.Template objects
        id_to_template: Dict[str, Template] = {} # mapping pyjess.Template ID to my Template objects
        for template_tuple in template_tuples:
            if check_template(template_tuple, warn=warn): # returns True if the Template passed all checks or if warn is set to False
                pyjess_template = template_tuple[0].to_pyjess_template()
                pyjess_templates.append(pyjess_template)
                id_to_template[pyjess_template.id] = template_tuple[0]

        # iterate over all the query molecules passed
        for molecule_path in molecule_paths:
            _single_query_run(outdir=outdir, molecule_path=molecule_path, pyjess_templates=pyjess_templates, id_to_template=id_to_template, rmsd=rmsd, distance=distance, max_dynamic_distance=max_dynamic_distance, conservation_cutoff=conservation_cutoff)

        # TODO what should i retrun? a path to the output file? a bool if any hits were found?
        # return
        
if __name__ == "__main__":
    
    parser = argparse.ArgumentParser()
    parser.add_argument('-i', '--input', type=Path, help='File path of the query pdb file OR File path of a file listing multiple query pdb files seperated by linebreaks')
    parser.add_argument('-j', '--jess', nargs = 4, type=float, help='Jess space seperated parameters rmsd, distance, max_dynamic_distance, score_cutoff')
    parser.add_argument('-o', '--output', type=Path, help='Output Directory to which results should get written')
    parser.add_argument('-v', '--verbose', type=bool, default=False, help='If process information and time progress should be printed to the command line')
    parser.add_argument('-w', '--warn', type=bool, default=True, help='If warings about bad template processing or suspicous and missing annotations should be raised')
    args = parser.parse_args()
    
    query_path = args.input
    jess_params = [i for i in args.jess]
    outdir = args.output
    warn = args.warn
    verbose = args.verbose

    # jess parameters
    # we use different parameters for different template residue numbers - higher number more generous parameters
    rmsd = jess_params[0] # in Angstrom, typcically set to 2
    distance = jess_params[1] # in Angstrom between 1.0 and 1.5 - lower is more strict. This changes with template size
    max_dynamic_distance = jess_params[2] # if equal to distance dynamic is off: this option is currenlty dysfunctional
    conservation_cutoff = jess_params[3] # Reads B-factor column in .pdb files: atoms below this cutoff will be disregarded. Could be pLDDT for example

    matcher_run(query_path=query_path, rmsd=rmsd, distance=distance, max_dynamic_distance=max_dynamic_distance, conservation_cutoff=conservation_cutoff, outdir=outdir, warn=warn, verbose=verbose)