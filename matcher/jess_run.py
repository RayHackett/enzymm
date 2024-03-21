"""
written by Ray
script six: jess_run.py
jess run script without filtering

Outsources to:
    Jess
    filter_output.py
    common_functions.py

Inputs:
    file containing a batch of pdb files and their respective absolute paths seperated by linebreaks supplied by -i flag
    Name pattern of Lists of Jess Templates by residue numbers from 3-8 in variable template_list
    Reads ./pdbs/<Query>_spectrum_js.pdb
    EClist_cofactors_forRH.csv .csv file mapping EC numbers to Cofactors. Got this from Neera
    cath-v4_3_0.alphafold-v2.2022-11-22.tsv .tsv file mapping CATH domains in Alphafold structures
    Jess Parameters rmsd, distance, max_dynamic_distance, score_cutoff
    
Outputs:
    ./jess results from jess seperated by residue number of the templates
    ./jess_logs log files from jess runs
    ./jess/<batch>all_res_info.json file with all info and annotations for all hits to the queries in the batch
    prints some useful info about what results are in which file
"""

import sys
import glob
import subprocess
import re
import argparse
import json
from pathlib import Path
from typing import Optional, Union, TextIO, Iterator
from functools import cached_property
from dataclasses import dataclass, field

import pyjess # cython port of Jess to python package

from .jess.filter_hits import filter_hits
from .template import Template, load_templates
from .common_files import chunks, ranked_argsort

d = Path(__file__).resolve().parent

def jess_call(molecule: TextIO, templates: List[Template], rmsd: float, distance_cutoff: float, max_dynamic_distance: float): #TODO this should be plddt cutoff not max dyn dist
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
    query
    ** inherit properties from Hit
    ## make sure you have
    # Template.ec
    # Template.cath
    # Template.size
    # Template.orientation_vector
    # Template.multimeric
    # Template.relative_order
    # Template.mcsa_id
    # Template.cluster
    # for res in Template.residues:
    #     res.orientation_vector
    # Hit.rmsd
    # Hit.log-evalue

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

    angle_mean = angle_mean([res.orientation_vector for res in Template.residues] , [calc_residue_orientation(atom_triplet) for atom_triplet in residues])

    # # Lookup
    # Hit.uniprot_id
    # Hit.pdb_id
    # Hit.ec
    # Hit.cath
    # Hit.interpro # TODO
    # Hit.is_catalytic

    # cofactors # TODO associated via query or template ec



    
def (filtered_output):
        
    mydict = {}
        ########################### Query name parsing #####################################################################
        query_name_splitter = re.split(r'\-|_|\.', Path(Remark_line.split()[1]).name) # second item is filepath of the query; we get the filename
        
        Uniprot_ID = ''
        for index, part in enumerate(query_name_splitter):
            # Check if Query is a Uniprot Accession number
            # see https://www.uniprot.org/help/accession_numbers
            match = re.search('[OPQ][0-9][A-Z0-9]{3}[0-9]|[A-NR-Z][0-9]([A-Z][A-Z0-9]{2}[0-9]){1,2}', part)
            if match:
                Uniprot_ID = match.group()
                Uniprot_index = index
        
        Query_is_af_CATH = False
        Query_is_PDBchain = False
        if Uniprot_ID and Uniprot_index != 0:
            # if the query is from Alphafold it has the structure [af|AF][\-|_][Uniprot_id]_xxxxx
            if query_name_splitter[Uniprot_index-1].upper() == 'AF': # if the bit before the Uniprot id is AF
                    # CATH Alphafold domains have the structure af_Uniprotid_start_finish
                    if len(query_name_splitter) == 4 and query_name_splitter[2].isdigit() and query_name_splitter[3].isdigit():
                        Query_is_af_CATH = True
                        Query = Path(Remark_line.split()[1]).name # filename of the query file
        
        # detecting pdb codes is currently a bit difficult as the 4 character ids dont have a particular pattern.
        # I am pretty sure though that the first character is always 1-9. Zero is not used.
        # In future this will not work and entry ids should be called pdb_00001abc_xyz_v1-2.cif
        # I assume then it is a pdb identifier with a chain ID
        elif re.search('[1-9]', query_name_splitter[0][0]): # check if the first character is in 1-9
            Query_is_PDBchain = True
            Query = query_name_splitter[0]
        else:
            Query = Path(Remark_line.split()[1]).name



        # if template and match are derived from the same pdb structure, do something
        if Query_is_PDBchain:
            if Template.pdb_id.upper() == Query.pdb_id.upper():
                continue # this skips this iteration and excludes the match
            
        ############################# Below Query Annoations relying on extral data sources ###################
        
        Retrieve:
        Query.uniprot_id
        Query.pdb_id
        Query.ec: set()
        Query.cath: set()

        if Query.pdb_chain:
            #update via sifts

        if Query.uniprot_id
            Query_EC.update(get_ec(Uniprot_ID)) # get EC from Uniprot
            
        # caths associated which each alphafold structure are loaded into dataframe cath_df
        all_query_caths = [] # this list will contain dicts with just one key value pair each
        if Query_is_af_CATH:
            query_df = cath_df.loc[cath_df['domain_ID'] == Query]
            cath_id = query_df['sfam_id'].to_list() # get the cath id
            residue_pos = query_df['domain_ID'].to_list() # this column contains start-end residue positions
            for i in range(len(cath_id)):
                all_query_caths.append({cath_id[i]: residue_pos[i].split('_')[2:]}) # value is list with [start, end]
        elif Uniprot_ID:
            query_df = cath_df.loc[cath_df['Uniprot_ID'] == Uniprot_ID]
            cath_id = query_df['sfam_id'].to_list() # get the cath id
            residue_pos = query_df['domain_ID'].to_list() # this column contains start-end residue positions
            for i in range(len(cath_id)):
                all_query_caths.append({cath_id[i]: residue_pos[i].split('_')[2:]}) # value is list with [start, end]
        

        """
        # add a list of cofactors associated with each EC number
        cofactors = set()
        if Template_EC in cofactor_dict:
            cofactors.update(cofactor_dict[Template_EC])
        for i in Query_EC:
            if i in cofactor_dict:
                cofactors.update(cofactor_dict[i])
        
        # check if it is catalytic by looking for EC or GO annotations
        is_catalytic = False
        if Uniprot_ID
            # Is it an enzyme or nonenzyme?
            is_catalytic = catalytic_checker(Uniprot_ID) # check if the entry has the catalytic GO annotation
        
        # Interpro function has its own checker to correctly handle UniProt or PDB accessions.
        if Query_is_PDBchain:
            domains, fams, superfams = get_interpro(Query)
        elif Uniprot_ID:
            domains, fams, superfams = get_interpro(Uniprot_ID)
        else:
            domains = []
            fams = []
            superfams = []
            
        Query_InterPro = {'domains': domains, 'families': fams, 'superfamilies': superfams}
        
        # add annotation if there is an active site or binding site annotation in Uniprot
        active_resids = []
        binding_resids = []
        if Uniprot_ID:
            active_resids, binding_resids = annotate_sites_from_Uniprot(Uniprot_ID)
            
        """
            
        #####################################################################################################


        
    # only after this dictionary is complete could we add the completeness tag
    # we need a list of all templates for each M-CSA Query pair
    # thus Group Hits by Query
    # For each query check if all Templates assigned to the same cluster targeted that structure
    # TODO report statistics on this: This percentage of queries had a complete active site as reported by the completeness tag
    # Filter this by template clusters with >1 member of course or report seperately by the number of clustermembers
    # or say like: This template cluster was always complete while this template cluster was only complete X times out of Y Queries matched to one member

    # note that completeness is assigned on the cluster level. first digit is cluster number
    # old jess templates had only 1 cluster anyway
    # second number is the current cluster member number
    # get the total number templates in the group which is the third cluster digit

    # # For each query group
    #     # extract list of all templates which hit that query
    #     templates = set()
    #     for match in match_list:
    #         templates.add(Template)

    #     # group the templates by M-csa id

    #     # For each template group, check the third digit with is the total number of templates in that template cluster
    #         total_clustermembers = int(Template.cluster.split('.')[-1])

    #         check if all the cluster members with all 2nd digit identiers up to and including total_clustermembers are present in the group,
    #         if so, assign the completeness tag:
    #             Hit.complete = True
return mydict
            
        
def main(start_file: Path,  rmsd: float, distance: float, max_dynamic_distance: float, score_cutoff: float, outdir: Path, search_incrementally: bool):
    
    def read_line_by_line(file: Path) -> Iterator[Path]:
        try:
            for path in open(start_file):
                yield path
        except FileNotFoundError:
            print(f"Input file not found: {start_file}")

    for path in read_line_by_line(start_file):
        if path.exists():
            # check if its a .pdb or .ent file
            if path.suffix.lower() not in {".pdb", ".ent"}:
                raise ValueError(f'File {path} did not have expected .pdb or .ent extensions')
        else:
            raise FileNotFoundError(f'Molecule file {path} does not exist')

    if not outdir.is_dir():
        raise FileNotFoundError(f"'{outdir}' is not an existing directory.")

    ############## Loading globally used dictionaries ##################
    
    global Uniprot_json_data
    Uniprot_json_data = {} # Initalize globally for use in the get_complete_Uniprot() function
    
    global Uniprot_EC_dict
    Uniprot_EC_dict = {} # Initalize globally for use in the get_ec() function
    
    global InterPro_dict
    InterPro_dict = {} # Initialize globally for use in the get_interpro() function
    
    global pdb_to_cath_dict
    pdb_to_cath_dict = {} # Initialize globally for use in the get_cath_from_interpro() function
    
    global Uniprot_catalytic_GO_dict
    Uniprot_catalytic_GO_dict = {} # Initialize globally for use in the is_catalytic() function
    
    global Uniprot_site_annotation_dict
    Uniprot_site_annotation_dict = {} # Initialize globally for use in the annotate_sites_from_Uniprot() function
    
    ############# Loading Files with annotation data ######################
    
    # # I got this csv from Neera originally. Data is from 2020
    # cofactor_df = pd.read_csv(Path(d, '../Downloads/EClist_cofactors_forRH.csv'))
    # global cofactor_dict
    # cofactor_dict = {}
    # for index, row in cofactor_df.iterrows():
    #     if row['EC'] not in cofactor_dict:
    #         cofactor_dict[row['EC']] = []
    #     cofactor_dict[row['EC']].append(row['Cof_ID'])
    
    # # load the tsv file from cath which maps cath domains to alphafold stuctures
    # # requires pandas
    # global cath_df
    # # ftp ftp://orengoftp.biochem.ucl.ac.uk/
    # # /alphafold/cath-v4.3.0-model-organisms/cath-v4_3_0.alphafold-v2.2022-11-22.tsv
    # cath_df = pd.read_csv(Path(d, '../Downloads/cath-v4_3_0.alphafold-v2.2022-11-22.tsv'), sep='\t')
    # domains = cath_df['domain_ID'].to_list()
    # Uniprot_IDs = []
    # for i in domains:
    #     Uniprot_IDs.append(i.split('_')[1])
    # cath_df['Uniprot_ID'] = Uniprot_IDs # add a column with just the Uniprot id
    
    # # This maps PDBchains to Uniprot IDs, CATHs and ECs
    # global pdb_sifts_df
    # pdb_sifts_df = pd.read_csv(Path(d, '../Downloads/pdb_sifts.csv'))
    # # PDBchain column has PDB in lowercase!
        
    ################# Running Jess ###################################
    print('jess parameters: ', rmsd, distance, max_dynamic_distance, score_cutoff)
    print('writing results to {}'.str(Path(outdir).absolute()))

    all_res_dict = {}
    
    # we iterate over templates starting with the ones with the largest ones
    # we start the search with all the input queries
    # by changing the working_file during the loop we can exclude queries from searches with smaller templates
    working_file = start_file # start_file is the original input of query structures to search
    scanned_molecule_paths = set()
    for template_res_num in [8, 7, 6, 5, 4, 3]:
        print('Now on template res num', template_res_num)
        # load the templates with a given residue number, yields tuple with (Template,Path)
        templates = list(load_templates(files(__package__).joinpath('jess_templates_20230210', '{}_residues'.format(template_res_num)))[0])
        
        # # Optional Template checking here
        # for template in templates:
        #     check_template(template) #if specific warnings are found, drop this template

        scanned_molecules = set()
        # Pass all the molecule paths as a list to the jess_call function
        for molecule_path in read_line_by_line(start_file):

            # if false, if hits for a structure were found with a large template, no not continue searching with smaller templates
            if complete_search == False:
                if molecule_path not in scanned_molecule_paths:
                    with open(molecule_path, 'r') as molecule:
                        # one could filter templates here too by any property or attribute of the template class in template.py
                        best_hit = jess_call(molecule, templates, rmsd_threshold=rmsd, distance_cutoff=distance, max_dynamic_distance=max_dynamic_distance)
                        if best_hit:
                            scanned_molecules.add(molecule_path)
            else:
                # search with smaller templates too even if searches with larger templates returned hits
                best_hit = jess_call(molecule, templates, rmsd_threshold=rmsd, distance_cutoff=distance, max_dynamic_distance=max_dynamic_distance)

        scanned_molecule_paths.extend(scanned_molecules)

        # mydict = parse_jess(filtered_output)
        # # now the results from mydict get added to the all_res_dict
        # # all the queries for which no matches were found get added to the new_working_list
        # new_working_list = []
        # for key, value_keys in mydict.items():
        #     if key not in all_res_dict:
        #         all_res_dict[key] = {}
        #     for sub_key, match_list in value_keys.items():
                
        #         # if a match was not found, add the query to the new_working_list
        #         for query in old_input_queries:
        #             if sub_key not in query:
        #                 new_working_list.append(query)
                
        #         if sub_key not in all_res_dict[key]:
        #             all_res_dict[key][sub_key] = []
        #         for match in match_list:
        #             all_res_dict[key][sub_key].append(match)
        
        # # queries for which no results with a given template size were found
        # new_working_file = Path(jess_path, 'working_queries_after_' + str(template_res_num) + 'res.txt')
        # with open(new_working_file, 'w') as f:
        #     for i in new_working_list:
        #         f.write(i + '\n')
        
        # with open(Path(jess_path,Path(filtered_outfile).name.split('.')[0] + '_info.json'), 'w') as f:
        #     json.dump(mydict, f)
        
    # with open(Path(jess_results_path,Path(working_file).name.split('.')[0] + 'all_res_info.json'), 'w') as f:
    #     json.dump(all_res_dict, f)
    
        
if __name__ == "__main__":
    
    parser = argparse.ArgumentParser()
    parser.add_argument('-i', '--input', type=Path, help='file containing a list of pdb files seperated by linebreaks.')
    parser.add_argument('-j', '--jess', nargs = '+', help='Jess space seperated parameters rmsd, distance, max_dynamic_distance, score_cutoff, optional flags as a string')
    parser.add_argument('-o', '--output', type=Path, help='Output Directory to which results should get written', default=Path.cwd())
    parser.add_argument('-c', '--complete', tuype=bool, help='If True continue search with smaller templates even after hits with larger template have been found', default=True)
    # TODO add complete_search as an argument too
    args = parser.parse_args()
    
    start_file = args.input
    outdir = args.output
    jess_params = [i for i in args.jess]

    # jess parameters
    # we use different parameters for different template residue numbers - higher number more generous parameters
    rmsd = jess_params[0] # in Angstrom, typcically set to 2
    distance = jess_params[1] # in Angstrom between 1.0 and 1.5 - lower is more strict. This changes with template size
    max_dynamic_distance = jess_params[2] # if equal to distance dynamic is off: this option is currenlty dysfunctional
    score_cutoff = jess_params[3] # Reads B-factor column in .pdb files: atoms below this cutoff will be disregarded. Could be pLDDT for example
    
    if len(jess_params) != 4:
        sys.exit('Wrong number of Jess Parameters') # TODO what kind of error should I throw?

    main(start_file, rmsd, distance, max_dynamic_distance, score_cutoff, outdir, complete_search=True)