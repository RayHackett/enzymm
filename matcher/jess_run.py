"""
Inputs:
    file containing a batch of pdb files and their respective absolute paths seperated by linebreaks supplied by -i flag
    Jess Parameters rmsd, distance, max_dynamic_distance, score_cutoff
    
Outputs:
    ./jess results from jess seperated by residue number of the templates
    ./jess_logs log files from jess runs
    ./jess/<batch>all_res_info.json file with all info and annotations for all hits to the queries in the batch
    prints some useful info about what results are in which file
"""

import argparse
from pathlib import Path
from typing import List, Optional, Union, TextIO, Iterator
from functools import cached_property
from dataclasses import dataclass, field
import csv

import pyjess # cython port of Jess to python package

from .jess.filter_hits import filter_hits
from .template import Template, load_templates
from .utils import chunks, ranked_argsort

def jess_call(molecule: TextIO, templates: List[Template], rmsd: float, distance_cutoff: float, max_dynamic_distance: float) -> List[pyjess.Hit]: #TODO this should be plddt cutoff not max dyn dist
    # molecule is a pdb object
    # templates is a list of template objects
    jess = Jess(templates)
    # Create a Jess instance and use it to query a molecule (a PDB structure) against the stored templates:
    for molecule in molecules:
        query = jess.query(molecule=molecule, rmsd_threshold=2.0, distance_cutoff=3.0, max_dynamic_distance=3.0)

        if query: # if any hits were found
            # retain only the hit with the lowest e-value for each query
            best_hits = filter_hits(query)

            # #The hits are computed iteratively, and the different output statistics are computed on-the-fly when requested:
            # for hit in best_hits:
            #     print(hit.molecule.id, hit.template.name, hit.rmsd, hit.log_evalue)
            #     for atom in hit.atoms():
            #         print(atom.name, atom.x, atom.y, atom.z)

        return best_hits

@dataclass
class Match:
    """"Class for storing annotated Jess hits"""
    template
    query # This is the target structure
    # ** inherit properties from Hit
    ## make sure you have
    # Template.ec
    # Template.cath
    # Template.size
    # Template.true_size
    # Template.multimeric
    # Template.relative_order
    # Template.mcsa_id
    # Template.cluster
    # for res in Template.residues:
    #     res.orientation_vector
    # Hit.rmsd
    # Hit.log-evalue
    # Hit.complete # initialize as false but change it later

    @cached_property
    def multimeric(self) -> bool:
        """`bool`: Boolean if the atoms in the hit stem from multiple protein chains
        """
        return all(atom.chain == self.atoms[0].chain for atom in self.atoms)

    @cached_property
    def preserved_resid_order(self) -> bool:
        """`bool`: Boolean if the residues in the template and in the matched query structure have the same relative order.
        This is a good filtering parameter but excludes hits on examples of convergent evolution or circular permutations
        """
        if self.multimeric:
            return False
        else:
            # Now extract relative atom order in hit
            return ranked_argsort([atom.residue_number for atom in self.atoms]) == Template.relative_order

    @classmethod
    def calc_residue_orientation(cls, atoms: Tuple[pyjess.Atom, pyjess.Atom, pyjess.Atom]) -> Template.Vec3:
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

        vectup = vector_atom_type_dict[atoms[0].residue_name]
        # In residues with two identical atoms, the vector is calculated between the middle atom and the mid point between the identical pair
        if vectup[1] == 'mid':
            try:
                middle = next(Atom for Atom in atoms if Atom.name == vectup[0])
                side1, side2 = [atom for atom in atoms if atom != middle]
                midpoint = [(side1.x + side2.x) / 2, (side1.y + side2.y) / 2, (side1.z + side2.z) / 2]
                return Vec3(middle.x - midpoint[0], middle.y - midpoint[1], middle.z - midpoint[2])
            except StopIteration:
                raise ValueError(f"Failed to find middle atom for amino-acid {atoms[0].residue_name!r}") from None
                
        else: # from first atom to second atom
            try:
                first_atom = next(Atom for Atom in atoms if Atom.name == vectup[0])
            except StopIteration:
                raise ValueError(f'Failed to find first atom for amino-acid {atoms[0].residue_name!r}') from None
            try:
                second_atom = next(Atom for Atom in atoms if Atom.name == vectup[1])
            except StopIteration:
                raise ValueError(f'Failed to find second atom for amino-acid {atoms[0].residue_name!r}') from None
                
            return first_atom - second_atom

    # TODO write a list with matched residues
    # # Hit.atoms is a list of matched atoms with all info on residue numbers and residue chain ids and atom types, ideally this should conserve order!!!
    residues : List[Tuple[pyjess.Atom, pyjess.Atom, pyjess.Atom]] = []
    for atom_triplet in chunks(Hit.atoms, 3): # yield chunks of 3 atoms each
        # check if all three atoms belong to the same residue by adding a tuple of their residue defining properties to a set
        unique_residues = { (atom.residue_name, atom.chain_id, atom.residue_number) for atom in atom_triplet }
        if len(unique_residues) != 1:
            raise ValueError('Mixed up atom triplets. The atoms come from different residues!')
        residues.append(atom_triplet)
    
    @classmethod
    def angle_mean(temp_vec_list: List[Template.Vec3], match_vec_list: List[Template.Vec3]):
        if len(temp_vec_list) != match_vec_list:
            raise ValueError('Vector lists for Template and matching Query structure had different lengths.')

        def angle_between(v1: Template.Vec3, v2: Template.Vec3): # from https://stackoverflow.com/questions/2827393/angles-between-two-n-dimensional-vectors-in-python
            """ Returns the angle in radians between vectors 'v1' and 'v2':
            """
            v1_u = vector / np.linalg.norm(v1) # unit vector
            v2_u = vector / np.linalg.norm(v2) # unit vector
            return np.arccos(np.clip(np.dot(v1_u, v2_u), -1.0, 1.0))

        # now calculate the angle between the vector of the template and the query per residue
        angle_list = []
        for i in range(len(temp_vec_list)):
            angle_list.append(angle_between(match_vec_list[i], temp_vec_list[i]))

        return np.mean(angle_list)

    orientation = angle_mean([res.orientation_vector for res in Template.residues] , [calc_residue_orientation(atom_triplet) for atom_triplet in residues])

def matcher_run(molecule_path: Path, rmsd: float, distance: float, max_dynamic_distance: float, score_cutoff: float, outdir: Path, complete_search: bool):
    # TODO technically I want the user to be able to pass a filepath too so not sure about typing here

    if not outdir.is_dir():
        # TODO would it be better to create the directory instead and print a statement that a directory was created?
        raise FileNotFoundError(f"'{outdir}' is not an existing directory.")
    
    # TODO using Path(molecule_path) here just in case molecule_path is a string
    if Path(molecule_path).exists():
        # check if its a .pdb or .ent file
        if Path(molecule_path).suffix.lower() not in {".pdb", ".ent"}:
            raise ValueError(f'File {path} did not have expected .pdb or .ent extensions')

        ################# Running Jess ###################################
        # TODO set level of verbosity and then print if conditions apply
        print('jess parameters: ', rmsd, distance, max_dynamic_distance, score_cutoff)
        print('complete_search is set to: ', complete_search)
        print('writing results to {}'.str(Path(outdir).absolute()))
        
        # we iterate over templates starting with the ones with the largest ones
        all_hits = []
        for template_res_num in [8, 7, 6, 5, 4, 3]:
            print('Now on template res num', template_res_num)
            # load the templates with a given residue number, yields tuple with (Template,Path)
            templates = list(load_templates(files(__package__).joinpath('jess_templates_20230210', '{}_residues'.format(template_res_num)))[0])
            
            # # TODO Optional Template checking here
            # for template in templates:
            #     check_template(template) #if specific warnings are found, drop this template
            # one could filter templates here too by any property or attribute of the template class in template.py

            # if false, if hits for a structure were found with a large template, no not continue searching with smaller templates
            if complete_search == False:
                if not all_hits:
                    with open(molecule_path, 'r') as molecule:
                        best_hits = jess_call(molecule, templates, rmsd_threshold=rmsd, distance_cutoff=distance, max_dynamic_distance=max_dynamic_distance)
                    if best_hits:
                        all_hits.extend(best_hits)
            else: # search with smaller templates too even if searches with larger templates returned hits
                with open(molecule_path, 'r') as molecule:
                    best_hist = jess_call(molecule, templates, rmsd_threshold=rmsd, distance_cutoff=distance, max_dynamic_distance=max_dynamic_distance)
                if best_hits:
                    all_hits.extend(best_hits)

            # TODO what kind of method is this? is it even a function? 
            def _check_completeness():
                # only after all templates of a certain size have been scanned could we compute the completeness tag
                
                # Group hit objects by the pair (hit.template.m-csa, first digit of hit.template.cluster)
                def get_key(obj) -> Tuple[str, str]:
                    return obj.template.msca_id, Template.cluster.split('.')[0] # split on . and take the first digit which is the cluster number

                grouped_hits = [list(g) for _, g in itertools.groupby(sorted(all_hits, key=get_key), get_key)]

                for cluster_hits in grouped_hits:
                    # For each query check if all Templates assigned to the same cluster targeted that structure
                    #
                    # TODO report statistics on this: This percentage of queries had a complete active site as reported by the completeness tag
                    # Filter this by template clusters with >1 member of course or report seperately by the number of clustermembers
                    # or say like: This template cluster was always complete while this template cluster was only complete X times out of Y Queries matched to one member
                    #
                    # get the total number templates from the first hit in each cluster which is the third cluster digit
                    total_clustermembers = int(cluster_hits[0].Template.cluster.split('.')[2])
                    # check if all the cluster members with all 2nd digit identiers up to and including total_clustermembers are present in the group,
                    indexed_possible_cluster_members = list(range(total_clustermembers))
                    possible_cluster_members = [x+1 for x in indexed_possible_cluster_members]
                    
                    found_cluster_members = [int(hit.Template.cluster.split('.')[1]) for hit in cluster_hits] # second number is the current cluster member number

                    if found_cluster_members == possible_cluster_members:
                        for hit in cluster_hits:
                            hit.complete = True

        # TODO fix writing for object attributes/properties that are lists
        # TODO should results dir be written by the pipeline?
        p = Path(outdir, "results/")
        p.mkdir(parents=True, exist_ok=True)
        # TODO somehow we need to give unique filenames so stuff doesnt get overwritten
        with open(Path(p, 'query_matches.tsv'), 'w', newline='', encoding ="utf-8") as tsvfile:
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
            for hit in all_hits: # one line per match
                writer.writerow([
                                hit.query.id,
                                hit.template.id,
                                hit.template.size,
                                hit.template.true_size,
                                hit.template.mcsa_id,
                                hit.template.uniprot_id,
                                hit.template.pdb_id,
                                hit.template.ec,
                                hit.template.cath,
                                hit.template.multimeric,
                                hit.query.uniprot_id,
                                hit.query.pdb_id,
                                hit.query.ec,
                                hit.query.cath,
                                hit.query.multimeric,
                                hit.rmsd,
                                hit.log_evalue,
                                hit.orientation,
                                hit.preserved_order,
                                hit.completeness,
                                hit.matched_residues])
            
        # TODO what should i retrun? a path to the output file? a bool if any hits were found?
        # return

    else:
        raise FileNotFoundError(f'Molecule file {path} does not exist')
        
if __name__ == "__main__":
    
    parser = argparse.ArgumentParser()
    parser.add_argument('-i', '--input', type=Path, help='File path of the query pdb file.')
    parser.add_argument('-j', '--jess', nargs = '+', help='Jess space seperated parameters rmsd, distance, max_dynamic_distance, score_cutoff, optional flags as a string')
    parser.add_argument('-o', '--output', type=Path, help='Output Directory to which results should get written', default=Path.cwd())
    parser.add_argument('-c', '--complete', tuype=bool, help='If True continue search with smaller templates even after hits with larger template have been found', default=True)
    args = parser.parse_args()
    
    molecule_path = args.input
    jess_params = [i for i in args.jess]
    outdir = args.output
    complete_search = args.complete

    # jess parameters
    # we use different parameters for different template residue numbers - higher number more generous parameters
    rmsd = jess_params[0] # in Angstrom, typcically set to 2
    distance = jess_params[1] # in Angstrom between 1.0 and 1.5 - lower is more strict. This changes with template size
    max_dynamic_distance = jess_params[2] # if equal to distance dynamic is off: this option is currenlty dysfunctional
    score_cutoff = jess_params[3] # Reads B-factor column in .pdb files: atoms below this cutoff will be disregarded. Could be pLDDT for example
    
    if len(jess_params) != 4:
        raise ValueError('Wrong number of Jess Parameters')

    matcher_run(molecule_path=molecule_path, rmsd=rmsd, distance=distance, max_dynamic_distance=max_dynamic_distance, score_cutoff=score_cutoff, outdir=outdir, complete_search=complete_search)