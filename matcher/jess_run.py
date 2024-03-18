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

d = Path(__file__).resolve().parent

def jess_call(molecule: TextIO, templates: List[Template], rmsd: float, distance_cutoff: float, max_dynamic_distance: float): #TODO this should be plddt cutoff not max dyn dist
    # molecule is a pdb object
    # templates is a list of template objects
    # Create a Jess instance and use it to query a molecule (a PDB structure) against the stored templates:
    for molecule in molecules:
        query = jess.query(molecule, rmsd_threshold=2.0, distance_cutoff=3.0, max_dynamic_distance=3.0)

        if query: # if any hits were found
            # retain only the hit with the lowest e-value for each query
            best_hits = filter_hits(query)

            #The hits are computed iteratively, and the different output statistics are computed on-the-fly when requested:
            for hit in best_hits:
                print(hit.molecule.id, hit.template.name, hit.rmsd, hit.log_evalue)
                for atom in hit.atoms():
                    print(atom.name, atom.x, atom.y, atom.z)

        return best_hits

@dataclass
def Match:
    """"`int`: Class for storing annotated Jess hits"""
    template
    query
    ** inherit properties from hit

    # calculate
    avg_orientation
    preserved_residue_order

    # Lookup
    uniprot_id
    query_ec
    query_cath
    query_interpro
    query_is_catalytic

    #Lookup
    template_cath
    template_interpro

    cofactors # associated via query or template ec

    **all other template properties and attributes



    
def (filtered_output):
        
    mydict = {}
    for entry in range(len(split_entries)): # iterating over all matches; entry is now an index
        Remark_line = split_entries[entry][0] # first line in a match is the Remark line
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
                 
        # reading info from the remark line of the match
        MCSA_entry = int(re.split(r'_|\.', Path(split_entries[entry][0].split()[3]).name)[1])
        Template_EC = {ec_dict[int(MCSA_entry)]} # from the M-CSA
        Template = Remark_line.split()[3]
        rmsd = float(Remark_line.split()[2])
        D_val = float(Remark_line.split()[5])
        E_val = float(Remark_line.split()[7])
        Result = split_entries[entry] # the coordinate info
        
            
            
        # add additional EC and CATH annotations via sifts directly from the pdb from which the template was made
        Template_PDB = Path(Template).name.split('.')[2].split('_')[0]
        for chain in chains:
            Template_PDBchain = Template_PDB + chain
            temp_uniprot, temp_cath, temp_ec = pdb_from_sifts(Template_PDBchain)
            if temp_ec:
                Template_EC.update(temp_ec)
            if temp_cath:
                Template_CATH.update(temp_cath)
            
        # get interpro annotations for the Template at the superfamily, family and domain levels
        Template_InterPro = MCSA_interpro_dict[MCSA_entry]
            
        # get list of all residues of the query which were matched
        matched_resids = []
        residue_numbers = []
        for line in Result[1:-1]: # excluding REMARK and ENDMDL
            # chain in [20:22] - jess accepts two letter chain ids!
            matched_resids.append(str(line[20:22]).strip() + str(int(line[22:26])))
            residue_numbers.append(int(line[22:26]))
        matched_resids = list(dict.fromkeys(matched_resids)) #preserves order
        residue_numbers = list(dict.fromkeys(residue_numbers)) #preserves order

        # check if the order of residues in the Query matches that of the Template
        # this is a good filtering parameter but excludes hits on examples of convergent evolution or circular permutations
        Preserved_resid_order = False
        query_res_rel_order = [sorted(residue_numbers).index(i) for i in residue_numbers]
        if query_res_rel_order == Template_res_order_dict[Template]:
            Preserved_resid_order = True
        
        # calculate the average angle between the vectors of corresponding matched residues
        angle_mean = calculate_residue_orientation(Template, Result)
        
        # if template and match are derived from the same pdb structure, do something
        if Query_is_PDBchain:
            if Template_PDB[:4].upper() == Query[:4].upper():
                continue # this skips this iteration and excludes the match
            
        ############################# Below Annoations relying on extral data sources ###################
        
        Query_CATH = set()
        Query_EC = set()
        if Query_is_PDBchain: # we look up data from a dataframe which contains mappings of pdbchains from sifts
            Uniprot_ID, CATH_NUMBER, EC_NUMBER = pdb_from_sifts(Query)
            if EC_NUMBER:
                Query_EC.update(EC_NUMBER)
            if CATH_NUMBER:
                Query_CATH.update(CATH_NUMBER)
                
        if Uniprot_ID:
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
        
        # We only want to include CATHs that include residues which were matched
        # M-CSA only includes the catalytic caths too
        for i in all_query_caths:
            cath = list(i.keys())[0]
            for res in matched_resids:
                resnum = re.findall(r'\d+', res)[0] # extract the resnumber from the string with the chain id
                # if any residue falls into the residue range annotaed with a cath, add that cath
                if int(resnum) in range(int(i[cath][0]), int(i[cath][1])+1):
                    Query_CATH.add(cath)
                    break # matching one residue within the domain is sufficient to add the annotation

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

        if MCSA_entry not in mydict:
            mydict[MCSA_entry] = {}
        if Query not in mydict[MCSA_entry]:
            mydict[MCSA_entry][Query] = []
        mydict[MCSA_entry][Query].append({'Query': Query,
                                          'Uniprot_ID': Uniprot_ID,
                                          'Template': Template,
                                          'Temp_res_num': Temp_res_num,
                                          'Multimeric': multimeric,
                                          'RMSD': rmsd,
                                          'D_val': D_val,
                                          'E_val': E_val,
                                          #'Catalytic': is_catalytic,
                                          'Template_EC': list(Template_EC),
                                          'Query_EC': list(Query_EC), # set to list so that exporting as .json works
                                          #'Cofactors': list(cofactors),
                                          'Template_CATH': list(Template_CATH),
                                          'Query_CATH': list(Query_CATH),
                                          'Template_InterPro': Template_InterPro,
                                          #'Query_InterPro': Query_InterPro,
                                          #'Active_site_resids': active_resids,
                                          #'Binding_site_resids': binding_resids,
                                          'drelSASA': drelSASA,
                                          'angle_mean': angle_mean,
                                          'Matched_resids': matched_resids,
                                          'Preserved_resid_order': Preserved_resid_order,
                                          'Full_Result': Result})
        
    # only after this dictionary is complete could we add the completeness tag
    # we need a list of all templates for each M-CSA Query pair
    for key, key_vals in mydict.items():
        for sub_key, match_list in key_vals.items():
            # extract list of all templates for a given query
            templates = set()
            for match in match_list:
                templates.add(Path(match['Template']).name)
                
            for match in match_list:
                # note that completeness is assigned on the cluster level. first digit is cluster number
                # old jess templates had only 1 cluster anyway
                # get the total number templates in the group which is the third cluster digit
                total_group = int(Path(match['Template']).name.split('.')[1][-1])
                # generate a general name by replacing the 2nd digit with ?
                general_name = re.sub(r"(cluster_.)_\d_", r"\1_?_", Path(match['Template']).name)

                # generate all other possible template names by replacement of the 2nd digit
                possible_names = set()
                for i in range(1, total_group+1):
                    possible_names.add(general_name.replace('?', str(i)))
                # only if all the possible template names are in the extracted list for a given query
                if possible_names.issubset(templates):
                    match['complete'] = True
                else:
                    match['complete'] = False
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
    
    # global ec_dict
    # # m-csa entry map to EC number
    # # every m-csa entry only has one EC number
    # with open(Path(d, '../Downloads/MCSA_EC_mapping.json'), 'r') as f:
    #     ec_dict = json.load(f)
    # # json always converts keys to strings. we want int type
    # ec_dict = {int(k):str(v) for k,v in ec_dict.items()}
    
    # # I got this csv from Neera originally. Data is from 2020
    # cofactor_df = pd.read_csv(Path(d, '../Downloads/EClist_cofactors_forRH.csv'))
    # global cofactor_dict
    # cofactor_dict = {}
    # for index, row in cofactor_df.iterrows():
    #     if row['EC'] not in cofactor_dict:
    #         cofactor_dict[row['EC']] = []
    #     cofactor_dict[row['EC']].append(row['Cof_ID'])
    
    # global mcsa_cath
    # # m-csa entry map to CATH id
    # # entry have have multiple CATH ids
    # with open(Path(d, '../Downloads/MCSA_CATH_mapping.json'), 'r') as f:
    #     mcsa_cath = json.load(f)
    # # json always converts keys to strings. we want int type
    # mcsa_cath = {int(k): v for k,v in mcsa_cath.items()}
    
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
    
    # global MCSA_interpro_dict
    # # dictonariy mapping M-CSA entries to Interpro Identifiers
    # # Interpro Acceccesions at the Domain, Family and Superfamily level
    # # are searched for the reference sequences of each M-CSA entry.
    # # Note that an M-CSA entry may have multiple reference sequences
    # with open(Path(d, '../Downloads/MCSA_interpro_dict.json'), 'r') as f:
    #     MCSA_interpro_dict = json.load(f)
    # # json always converts keys to strings. we want int type
    # MCSA_interpro_dict = {int(k): v for k,v in MCSA_interpro_dict.items()}(
        
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
        # load the templates with a given residue number
        templates = list(load_templates(files(__package__).joinpath('jess_templates_20230210', '{}_residues'.format(template_res_num))))

        already_scanned = set()
        # Pass all the molecule paths as a list to the jess_call function
        for molecule_path in read_line_by_line(start_file):

            # if false, if hits for a structure were found with a large template, no not continue searching with smaller templates
            if search_incrementally == False:
                if molecule_path not in scanned_molecule_paths:
                    with open(molecule_path, 'r') as molecule:
                        # one could filter templates here too by any property or attribute of the template class in template.py
                        best_hit = jess_call(molecule, templates, rmsd_threshold=rmsd, distance_cutoff=distance, max_dynamic_distance=max_dynamic_distance)
                        if best_hit:
                            already_scanned.add(molecule_path)
            else:
                # search with smaller templates too even if searches with larger templates returned hits
                best_hit = jess_call(molecule, templates, rmsd_threshold=rmsd, distance_cutoff=distance, max_dynamic_distance=max_dynamic_distance)

        scanned_molecule_paths.extend(already_scanned)

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
    # TODO add search_incrementally as an argument too
    args = parser.parse_args()
    
    start_file = args.input
    jess_params = [i for i in args.jess]

    # jess parameters
    # we use different parameters for different template residue numbers - higher number more generous parameters
    rmsd = jess_params[0] # in Angstrom, typcically set to 2
    distance = jess_params[1] # in Angstrom between 1.0 and 1.5 - lower is more strict. This changes with template size
    max_dynamic_distance = jess_params[2] # if equal to distance dynamic is off: this option is currenlty dysfunctional
    score_cutoff = jess_params[3] # Reads B-factor column in .pdb files: atoms below this cutoff will be disregarded. Could be pLDDT for example
    
    if len(jess_params) != 4:
        sys.exit('Wrong number of Jess Parameters') # TODO what kind of error should I throw?

    outdir = args.output
    main(start_file, rmsd, distance, max_dynamic_distance, score_cutoff, outdir, search_incrementally=True)