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
import subprocess
import requests
import re
import argparse
import json
import pandas as pd
import numpy as np
from pathlib import Path

import filter_output #yannis filtering script
from ..common_functions import json_extract
from ..common_functions import request_url
from ..common_functions import chunks
import calc_residue_orientation_jess as res_orientation

d = Path(__file__).resolve().parent

def jess_call(template_res_num, pdb_file, template_file, jess_params):
    out_name = str(jess_path) + '/' + Path(pdb_file).name.split('.')[0] + '_?res.pdb'.replace('?',str(template_res_num))
    
    # jess parameters
    # ultimately different params for different template residue numbers - higher number more generous params
    rmsd = jess_params[0] # in Angstrom, set to 2
    distance = jess_params[1] # in Angstrom between 1.0 and 1.5 - lower is more strikt
    max_dynamic_distance = jess_params[2] # if equal to distance dynamic is off: this option is dysfunctional
    score_cutoff = jess_params[3] # Reads B-factor column: atoms below this cutoff will be disregarded
    
    # the first four params are rmsd, distance, max_dynamic_distance and score_cutoff.
    # The (optional) fifth is a string controlling output format
    if len(jess_params) == 5:
        optional_flags = jess_params[4] # optional flags as a string without spaces
        script = ' '.join(['jess', str(template_file), str(pdb_file), str(rmsd), str(distance), str(max_dynamic_distance), str(score_cutoff), str(optional_flags), '>', str(out_name)])
    elif len(jess_params) == 4:
        script = ' '.join(['jess', str(template_file), str(pdb_file), str(rmsd), str(distance), str(max_dynamic_distance), str(score_cutoff), '>', str(out_name)])
    else:
        sys.exit('Wrong number of Jess Parameters')
    
    subprocess.check_call(script, shell=True)
    
    return pdb_file, out_name

def jess_filter(input_file):    # filter with yannis script
    with open(input_file, 'r') as f:
        data = f.readlines()
    filtered_output = filter_output.main(data)
    
    # gets written to same dir as input_file
    filtered_outfile = input_file.replace('.pdb', '_filtered.pdb')
    with open(filtered_outfile, 'w') as f:
        for i in filtered_output:
            f.write(str(i) + '\n')
        
    return filtered_outfile

def get_interpro(Identifier):
    
    if Identifier in InterPro_dict:
        interpro_domains = InterPro_dict[Identifier][0]
        interpro_families = InterPro_dict[Identifier][1]
        interpro_superfam = InterPro_dict[Identifier][2]
        
    else:
        checker = re.search('[OPQ][0-9][A-Z0-9]{3}[0-9]|[A-NR-Z][0-9]([A-Z][A-Z0-9]{2}[0-9]){1,2}', Identifier) 
        if checker:
            Uniprot_ID = checker.group()
            url = 'https://www.ebi.ac.uk/interpro/api/entry/InterPro/protein/UniProt/' + Uniprot_ID
        else:
            PDB = Identifier[:4]
            url = 'https://www.ebi.ac.uk/interpro/api/entry/interpro/structure/pdb/' + PDB
        
        interpro_domains = set()
        interpro_families = set()
        interpro_superfam = set()
        r = request_url(url, [200, 204])
        if r:
            data = r.json()
            # by scraping only the stuff associated with certain amino acid positions this would be more accurate
            # look not only at the metadata stuff. but that gets a bit more complicated
            for i in data['results']:   
                if i['metadata']['type'] == 'domain':
                    interpro_domains.add(i['metadata']['accession'])
                elif i['metadata']['type'] == 'family':
                    interpro_families.add(i['metadata']['accession'])
                elif i['metadata']['type'] == 'homologous_superfamily':
                    interpro_superfam.add(i['metadata']['accession'])
        else:
            print('No InterPro Response for:')
            print(url)
            
        InterPro_dict[Identifier]=[list(interpro_domains), list(interpro_families), list(interpro_superfam)]
    
    return list(interpro_domains), list(interpro_families), list(interpro_superfam)

def get_cath_from_interpro(PDBchain):
    all_query_cath = []
    if PDBchain in pdb_to_cath_dict:
        all_query_cath = pdb_to_cath_dict[PDBchain]
    else:
        # get cath annotations from Interpro. this is easier than fenangling the CATH api
        url = 'https://www.ebi.ac.uk/interpro/api/entry/cathgene3d/structure/pdb/' + PDBchain[:4]
        r = request_url(url, [200, 204])
        if r:
            all_query_cath = []
            data = r.json()
            for i in data['results']:
                for struc in i['structures']:
                    if struc['chain'] == PDBchain[4:]: # find the right chain
                        cath_id = i['metadata']['accession'].split(':')[1] # get the accession
                        for loc in struc['entry_protein_locations']:
                            for frag in loc['fragments']: # and then the location of the domain numbered by PDB res number
                                start = frag['start']
                                end = frag['end']
                                all_query_cath.append({cath_id: [start, end]})
        else:
            print(url)          
        pdb_to_cath_dict[PDBchain] = all_query_cath
    return all_query_cath

def pdb_from_sifts(PDBchain):
    # using sifts mapping
    PDBchain = PDBchain[:4].lower() + PDBchain[4:] # pdb_sifts_df has PDBchain in lowercase; chain is case sensitive
    if pdb_sifts_df['PDBchain'].str.contains(PDBchain).any(): # float('nan') is not equal to float('nan')
        CATH_set = set(pdb_sifts_df[pdb_sifts_df['PDBchain']==PDBchain]['CATH_NUMBER'].to_list())
        CATH_NUMBER = {x for x in CATH_set if x==x}
                
        EC_set = set(pdb_sifts_df[pdb_sifts_df['PDBchain']==PDBchain]['EC_NUMBER'].to_list())
        EC_NUMBER = {x for x in EC_set if x==x} # nan is float('nan') equal to float('nan')
        
        Uniprot_ID = pdb_sifts_df[pdb_sifts_df['PDBchain']==PDBchain]['UNIPROT_ID'].iloc[0]
        if pd.isna(Uniprot_ID):
            Uniprot_ID = ''
        return Uniprot_ID, CATH_NUMBER, EC_NUMBER
    else:
        return '' , '' , ''

def catalytic_checker(Uniprot_ID):
    if Uniprot_ID in Uniprot_catalytic_GO_dict:
        has_catalytic_go = Uniprot_catalytic_GO_dict[Uniprot_ID]
        
    else:
        has_catalytic_go = False
        try:
            url = 'https://rest.uniprot.org/uniprotkb/stream?format=list&query=%28%28accession%3A{}%29%20AND%20%28go%3A0003824%29%29'.format(Uniprot_ID)
            r = request_url(url, [200])
            if r.text:
                has_catalytic_go = True
            Uniprot_catalytic_GO_dict[Uniprot_ID] = has_catalytic_go
        except:
            pass
        
    return has_catalytic_go

def get_complete_Uniprot(Uniprot_ID):
    if Uniprot_ID in Uniprot_json_data:
        data = Uniprot_json_data[Uniprot_ID]
    else:
        url = 'https://rest.uniprot.org/uniprotkb/' + Uniprot_ID + '.json'
        r = request_url(url, [200])
        data = r.json()
        Uniprot_json_data[Uniprot_ID] = data
            
    return data

def get_ec(Uniprot_ID):
    if Uniprot_ID in Uniprot_EC_dict:
        ec_num = Uniprot_EC_dict[Uniprot_ID]
        
    else:
        data = get_complete_Uniprot(Uniprot_ID)
        ec_num = []
        ec_num.extend(json_extract(data, 'ecNumber'))
        for i in json_extract(data, 'ecNumbers'):
            for j in i:
                ec_num.append(j['value'])
        ec_num = list(set(ec_num))
                
        Uniprot_EC_dict[Uniprot_ID] = ec_num
        if not ec_num:
            print('No EC found for Uniprot Entry ' + Uniprot_ID)
    
    return ec_num
            
def annotate_sites_from_Uniprot(Uniprot_ID):
    if Uniprot_ID in Uniprot_site_annotation_dict:
        active_resids = Uniprot_site_annotation_dict[Uniprot_ID][0]
        binding_resids = Uniprot_site_annotation_dict[Uniprot_ID][1]
    else:
        data = get_complete_Uniprot(Uniprot_ID)
        active_resids = []
        binding_resids = []
        if 'features' in data:
            for feature in data['features']:
                if feature['type'] == 'Active site':
                    start = feature['location']['start']['value']
                    end = feature['location']['end']['value']
                    active_resids.extend(list(range(start,end+1)))

                elif feature['type'] == 'Binding site':
                    start = feature['location']['start']['value']
                    end = feature['location']['end']['value']
                    binding_resids.extend(list(range(start,end+1)))
                
        Uniprot_site_annotation_dict[Uniprot_ID] = [active_resids, binding_resids]
        
    return active_resids, binding_resids

def calculate_similarity_in_rsas(Template, Full_Result): # calculates solvent accessiblity differences
    # requires Template_res_dssp_dict which is calculated in the prerun step!
    # requires that solvent accessiblity values are saved in the B-factor column for the query structures!
    query_rsa_vals = []
    seen_resnums = set()
    for line in Full_Result[1:-1:3]: # reading every third line excluding the REMARK and ENDMDL line
        resname = line[17:20]
        chain = line[20:22].strip()
        resnum = int(line[22:26])
        if resnum not in seen_resnums: # only reading unique residue numbers as some residues occur twice
            query_rsa_vals.append(float(line[55:60])) # B-factor column!
        seen_resnums.add(resnum)
    
    # get the solvent accessibility values for the template from dictionary
    if Template_res_dssp_dict[Template]:
        template_rsa_vals = list(Template_res_dssp_dict[Template].values())
    else:
        print('lacking rsas for', Template)
        return 'NA'

    rsa_differences = []
    if len(template_rsa_vals) == len(query_rsa_vals): # need to have the same number of residues!
        for i in range(len(query_rsa_vals)):
            # DSSP puts 9.99 if the value could not be calculated
            if query_rsa_vals[i] != 9.99 and template_rsa_vals[i] != 9.99: 
                rsa_differences.append(abs(query_rsa_vals[i]-template_rsa_vals[i]))
                
    # only return a value if drelSASA for at least 3 residue pairs was calculated
    if len(rsa_differences) >= 3: 
        return np.sqrt(np.mean(np.square(np.array(rsa_differences))))
    else:
        return float('NaN')
    
def calculate_residue_orientation(Template, Full_Result):
    # Requires Template_vec_dict which is calculated in the prerun step!
    
    # for every residue in the match we calculate a vector and add it to a list
    residue_generator = chunks(Full_Result[1:-1], 3)
    query_vec_list = []
    for index, residue in enumerate(residue_generator):
        Template_residue = Template_vec_dict[Template][str(index)]
        vec_info = Template_residue['index_order']
        query_vec = res_orientation.get_query_vec(residue, vec_info)
        query_vec_list.append(query_vec)

    # calculate the angle equivalent of RMSD
    # note that this requires the same reference coordinate system for query and template
    template_vectors = [Template_vec_dict[Template][str(i)]['vector'] for i in range(len(query_vec_list))]
    template_residue_positions = Template_vec_dict[Template]['residue_positions']
    
    # angle_mean and angle_rms correlate the best with RMSD.
    # angle_mean is simpler to understand therefore we use this one
    mean_val = res_orientation.angle_mean(template_vectors, query_vec_list)
    
    #rms_val = res_orientation.angle_rms(template_vectors, query_vec_list)
    #angle_rmsd = res_orientation.angle_difference_test(template_vectors, query_vec_list)
    #weighted_angle_rmsd = res_orientation.weighted_angle_difference_test(template_vectors, query_vec_list, template_residue_positions)
    
    return mean_val
    
def parse_jess(filtered_file):
    # parse the output from jess
    with open(filtered_file, 'r') as f:
        data = f.read()
    entries = data.split('\n\n') # double linebreak seperarates matches
    split_entries = [] # list of matches
    for entry in entries:
        split_entry = entry.split('\n') # split_entry is a list of lines
        if all(split_entry):
            split_entries.append(split_entry) # list of lists
        
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
        
        # Temp_res_num is not necessarily equal to the resnumber given in path:
        # Residues matching ANY amino acid type are not counted as they are too general
        # Even if a Residue must match a particular AA, 6 atoms from a given residue may be selected
        # once for main chain and once for side chain
        # Therefore we only count unique template residues!
        Temp_res_num = Template_res_num_dict[Template]
            
        # check if the Template is split over multiple chains
        # this method will fail with chain names such as A1 which contain numeric characters
        res_list = Path(Template).name.split('.')[2].split('_')[1].split('-')
        chains = set(re.sub('[0-9]', '', i) for i in res_list) # get only the chain name by removing numeric characters
        if len(chains) > 1: # more than one chain means multimeric template
            multimeric = True
        else:
            multimeric = False
            
        # add Template CATH annotations from M-CSA
        if mcsa_cath[int(MCSA_entry)]:
            Template_CATH = set(mcsa_cath[int(MCSA_entry)])
        else:
            Template_CATH = set()
            
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
        for line in Result[1:-1]: # excluding REMARK and ENDMDL
            # chain in [20:22] - jess accepts two letter chain ids!
            matched_resids.append(str(line[20:22]).strip() + str(int(line[22:26])))
        matched_resids = list(dict.fromkeys(matched_resids)) #preserves order
        
        # check if the order of residues in the Query matches that of the Template
        # this is a good filtering parameter but excludes hits on examples of convergent evolution or circular permutations
        Preserved_resid_order = False
        query_res_rel_order = [sorted(matched_resids).index(i) for i in matched_resids]
        if query_res_rel_order == Template_res_order_dict[Template]:
            Preserved_resid_order = True
            
            
        # calculate difference in relative solvent accessivility between template and query structure
        drelSASA = calculate_similarity_in_rsas(Template, Result)
        
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
                general_name = re.sub("(cluster_.)_\d_", r"\1_?_", Path(match['Template']).name)

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

def merge_dicts(dica, dicb):
    dicc={}
    for key_a, value_a in dica.items():
        if key_a in dicb.keys():
            dicc[key_a] = dica[key_a] | dicb[key_a]
        else:
            dicc[key_a] = value_a
    for key_b, value_b in dicb.items():
        if key_b not in dicc.keys():
            dicc[key_b] = value_b
    return dicc
            
        
def main(start_file, jess_params):
    
    ########## Checking Paths and making folders ###############
    global cwd
    cwd = Path.cwd()
    
    # folder for the templates used
    global template_folder
    template_folder = Path(cwd, 'templates_used')
    template_folder.mkdir(exist_ok=True)
    
    # pdb path
    global pdb_path
    pdb_path = Path(cwd,'pdbs') 
    pdb_path.mkdir(exist_ok=True)
    
    # jess path
    global jess_path
    jess_path = Path(cwd,'jess') 
    jess_path.mkdir(exist_ok=True)
    
    # jess output path
    global jess_results_path
    jess_results_path = Path(cwd,'jess_results') 
    jess_results_path.mkdir(exist_ok=True)
    
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
    
    global ec_dict
    # m-csa entry map to EC number
    # every m-csa entry only has one EC number
    with open(Path(d, '../Downloads/MCSA_EC_mapping.json'), 'r') as f:
        ec_dict = json.load(f)
    # json always converts keys to strings. we want int type
    ec_dict = {int(k):str(v) for k,v in ec_dict.items()}
    
    # I got this csv from Neera originally. Data is from 2020
    cofactor_df = pd.read_csv(Path(d, '../Downloads/EClist_cofactors_forRH.csv'))
    global cofactor_dict
    cofactor_dict = {}
    for index, row in cofactor_df.iterrows():
        if row['EC'] not in cofactor_dict:
            cofactor_dict[row['EC']] = []
        cofactor_dict[row['EC']].append(row['Cof_ID'])
    
    global mcsa_cath
    # m-csa entry map to CATH id
    # entry have have multiple CATH ids
    with open(Path(d, '../Downloads/MCSA_CATH_mapping.json'), 'r') as f:
        mcsa_cath = json.load(f)
    # json always converts keys to strings. we want int type
    mcsa_cath = {int(k): v for k,v in mcsa_cath.items()}
    
    # load the tsv file from cath which maps cath domains to alphafold stuctures
    # requires pandas
    global cath_df
    # ftp ftp://orengoftp.biochem.ucl.ac.uk/
    # /alphafold/cath-v4.3.0-model-organisms/cath-v4_3_0.alphafold-v2.2022-11-22.tsv
    cath_df = pd.read_csv(Path(d, '../Downloads/cath-v4_3_0.alphafold-v2.2022-11-22.tsv'), sep='\t')
    domains = cath_df['domain_ID'].to_list()
    Uniprot_IDs = []
    for i in domains:
        Uniprot_IDs.append(i.split('_')[1])
    cath_df['Uniprot_ID'] = Uniprot_IDs # add a column with just the Uniprot id
    
    # This maps PDBchains to Uniprot IDs, CATHs and ECs
    global pdb_sifts_df
    pdb_sifts_df = pd.read_csv(Path(d, '../Downloads/pdb_sifts.csv'))
    # PDBchain column has PDB in lowercase!
    
    global MCSA_interpro_dict
    # dictonariy mapping M-CSA entries to Interpro Identifiers
    # Interpro Acceccesions at the Domain, Family and Superfamily level
    # are searched for the reference sequences of each M-CSA entry.
    # Note that an M-CSA entry may have multiple reference sequences
    with open(Path(d, '../Downloads/MCSA_interpro_dict.json'), 'r') as f:
        MCSA_interpro_dict = json.load(f)
    # json always converts keys to strings. we want int type
    MCSA_interpro_dict = {int(k): v for k,v in MCSA_interpro_dict.items()}
    
    # load dictionary mapping Template to number of residues
    global Template_res_num_dict
    with open(Path(d, '../Downloads/Template_res_num_dict.json'), 'r') as f:
        Template_res_num_dict = json.load(f)
    
    # load dictionary mapping Template to the relative order of residues therein
    global Template_res_order_dict
    with open(Path(d, '../Downloads/Template_res_order_dict.json'), 'r') as f:
        Template_res_order_dict = json.load(f)
        
    # load the dicionary with solvent accessible surface area for each residue in each template
    global Template_res_dssp_dict
    with open(Path(d, '../Downloads/Template_res_dssp_dict.json'), 'r') as f:
        Template_res_dssp_dict = json.load(f)
        
    # load the dictionary with the orientation vector for each residue in each template
    global Template_vec_dict
    with open(Path(d, '../Downloads/Template_vec_dict.json'), 'r') as f:
        Template_vec_dict = json.load(f)
        
    ################# Running Jess ###################################
    print('jess parameters: ', jess_params)
    # this is the final output dictionary
    all_res_dict = {}
    
    # we iterate over templates starting with the ones with the largest ones
    # we start the search with all the input queries
    # by changing the working_file during the loop we can exclude queries from searches with smaller templates
    working_file = start_file # start_file is the original input of query structures to search
    for template_res_num in [8, 7, 6, 5, 4, 3]:
        
        template_file = Path(template_folder, 'list_of_{}_residue_templates.list'.format(str(template_res_num)))
        print('Now on template res num', template_res_num)
        # prevous_file is the used input file (= the previous working file)
        # out_name is the name of the file to which the results of the jess search were written
        previous_file, out_name = jess_call(template_res_num, working_file, template_file, jess_params)
        print('jess input was', previous_file)
        print('jess writing to', out_name)
        
        # make a list of all the scanned queries
        with open(previous_file, 'r') as f:
            query_lines = f.readlines()
        old_input_queries = [i.strip() for i in query_lines]
        
        filtered_outfile = jess_filter(out_name)
        print('filtering was written to', filtered_outfile)
        mydict = parse_jess(filtered_outfile)
        # now the results from mydict get added to the all_res_dict
        # all the queries for which no matches were found get added to the new_working_list
        new_working_list = []
        for key, value_keys in mydict.items():
            if key not in all_res_dict:
                all_res_dict[key] = {}
            for sub_key, match_list in value_keys.items():
                
                # if a match was not found, add the query to the new_working_list
                for query in old_input_queries:
                    if sub_key not in query:
                        new_working_list.append(query)
                
                if sub_key not in all_res_dict[key]:
                    all_res_dict[key][sub_key] = []
                for match in match_list:
                    all_res_dict[key][sub_key].append(match)
        
        # queries for which no results with a given template size were found
        new_working_file = Path(jess_path, 'working_queries_after_' + str(template_res_num) + 'res.txt')
        with open(new_working_file, 'w') as f:
            for i in new_working_list:
                f.write(i + '\n')
                
        #######################################################
        # We can change the working_file to new_working_file 
        # to exclude already matched queries from future searches
        # simply uncomment the following line:
        # working_file = new_working_file
        
        # might need to change the way output files are aggregated at the end to make this work
        
        # OR we can just keep the start_file as our working file
        # if we dont want to exclude any queries from searches with smaller templates
        ########################################################
        
        with open(Path(jess_path,Path(filtered_outfile).name.split('.')[0] + '_info.json'), 'w') as f:
            json.dump(mydict, f)
        
    with open(Path(jess_results_path,Path(working_file).name.split('.')[0] + 'all_res_info.json'), 'w') as f:
        json.dump(all_res_dict, f)
        
if __name__ == "__main__":
    
    parser = argparse.ArgumentParser()
    parser.add_argument('-i', '--input', help='file containing a list of pdb files seperated by linebreaks.')
    parser.add_argument('-j', '--jess', nargs = '+', help='Jess space seperated parameters rmsd, distance, max_dynamic_distance, score_cutoff, optional flags as a string')
    args = parser.parse_args()
    
    start_file = args.input
    jess_params = [i for i in args.jess]
    main(start_file, jess_params)
    
