"""
written by Ray
jess_prerun.py

Generate and or check files necessary to run jess:

Outsources to:
    common_functions.py

Inputs:
    files with paths to each each template split by template residue number
    
Outputs:
    ./Downloads/MCSA_EC_mapping.json
    ./Downloads/MCSA_CATH_mapping.json
    ./Downloads/MCSA_interpro_dict.json
    ./Downloads/Template_res_num_dict.json
    ./Downloads/Resnum_per_Template_dict.json
    ./Downloads/Template_vec_dict.json
    
"""

import sys
from pathlib import Path
import json
import re
import numpy as np
import pandas as pd
import urllib.request
import gzip
import shutil

from .common_functions import request_url

def get_mcsa_ec():
    # m-csa entry map to EC number
    # every m-csa entry only has one EC number
    if not Path('./Downloads/MCSA_EC_mapping.json').is_file():
        ec_dict = {}
        for page in range(1,12):
            url = 'https://www.ebi.ac.uk/thornton-srv/m-csa/api/entries/?format=json&page='+str(page)
            general = request_url(url, [200])
            general_page = general.json()

            for i in range(len(general_page['results'])):
                mcsa_id = int(general_page['results'][i]['mcsa_id'])
                if mcsa_id not in ec_dict:
                    ec_dict[mcsa_id] = ''
                if len(general_page['results'][i]['all_ecs']) > 1:
                    print(mcsa_id)
                for ec_num in general_page['results'][i]['all_ecs']:
                    ec_dict[mcsa_id] = ec_num

            Path('./Downloads').mkdir(parents=True, exist_ok=True)
            with open(Path('./Downloads/MCSA_EC_mapping.json'), 'w') as f:
                json.dump(ec_dict, f)

def get_mcsa_cath():
    # m-csa entry map to CATH id
    # entry have have multiple CATH ids
    if not Path('./Downloads/MCSA_CATH_mapping.json').is_file():
        mcsa_cath = {}
        for page in range(1,12):
            url = 'https://www.ebi.ac.uk/thornton-srv/m-csa/api/entries/?format=json&page='+str(page)
            general = request_url(url, [200])
            general_page = general.json()
            
            for i in range(len(general_page['results'])):
                cath_ids = set()
                for resid in general_page['results'][i]['residues']:
                    for chain in resid['residue_chains']:
                        cath_ids.add(chain['domain_cath_id'])

                if '' in cath_ids:
                    cath_ids.remove('')
                cath_ids = list(cath_ids)
                mcsa_cath[general_page['results'][i]['mcsa_id']] = cath_ids
            
            Path('./Downloads').mkdir(parents=True, exist_ok=True)
            with open(Path('./Downloads/MCSA_CATH_mapping.json'), 'w') as f:
                json.dump(mcsa_cath, f)
                
def get_interpro(Identifier):
        
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
    
    return list(interpro_domains), list(interpro_families), list(interpro_superfam)

def get_mcsa_interpro():
    # dictonariy mapping M-CSA entries to Interpro Identifiers
    # Interpro Acceccesions at the Domain, Family and Superfamily level
    # are searched for the reference sequences of each M-CSA entry.
    # Note that an M-CSA entry may have multiple reference sequences
    if not Path('./Downloads/MCSA_interpro_dict.json').is_file():
        mcsa_references = {}
        for page in range(1,12):
            url = 'https://www.ebi.ac.uk/thornton-srv/m-csa/api/entries/?format=json&page='+str(page)
            general = request_url(url, [200])
            general_page = general.json()
            
            for i in range(len(general_page['results'])):
                mcsa_references[general_page['results'][i]['mcsa_id']] = general_page['results'][i]['reference_uniprot_id']
        
        MCSA_interpro_dict = {}

        for key, reference_seqs in mcsa_references.items():
            reference_seqs = reference_seqs.split(',')
            template_interpro_domains = set()
            template_interpro_fams = set()
            template_interpro_superfams = set()
            for seq in reference_seqs:
                seq = seq.replace(' ','')
                domains , fams, superfams = get_interpro(seq)
                template_interpro_domains.update(domains)
                template_interpro_fams.update(fams)
                template_interpro_superfams.update(superfams)

            MCSA_interpro_dict[key] = {'domains': list(template_interpro_domains), 'families': list(template_interpro_fams), 'superfamilies': list(template_interpro_superfams)}
        
        Path('./Downloads').mkdir(parents=True, exist_ok=True)
        with open(Path('./Downloads/MCSA_interpro_dict.json'), 'w') as f:
                json.dump(MCSA_interpro_dict, f)
                
def make_pdb_sifts_df():
    if not Path('./Downloads/pdb_sifts.csv').is_file():
        check_other_files()

        PDBchain_to_CATH_Uniprot = Path('./data/pdb_chain_cath_uniprot.csv')
        PDBchain_to_EC_Uniprot = Path('./data/pdb_chain_enzyme.csv')
        CATH_names_mapping = Path('./Downloads/cath-domain-list.txt')

        # load them all as df
        pdb_cath_uniprot_df = pd.read_csv(PDBchain_to_CATH_Uniprot, comment='#', names=['PDB', 'CHAIN', 'UNIPROT_ID', 'CATH_ID'])
        pdb_enzyme_uniprot_df = pd.read_csv(PDBchain_to_EC_Uniprot, comment='#', names=['PDB', 'CHAIN', 'UNIPROT_ID', 'EC_NUMBER'])
        cath_naming_df = pd.read_csv(CATH_names_mapping, comment='#', sep='\s+', header=None, dtype=str)

        # join the right columns to form the CATH NUMBER
        cath_naming_df['CATH_NUMBER'] = cath_naming_df[1] + '.' + cath_naming_df[2] + '.' + cath_naming_df[3] + '.' + cath_naming_df[4]
        # rename a column
        cath_naming_df.rename(columns = {0:'CATH_ID'}, inplace=True)

        # merge the three dataframes
        temp = pd.merge(pdb_cath_uniprot_df, pdb_enzyme_uniprot_df, how="outer", on=None) # outer merge on the first two
        pdb_sifts_df = temp.merge(cath_naming_df[['CATH_ID', 'CATH_NUMBER']], on='CATH_ID') # add in the CATH number by CATH_ID
        pdb_sifts_df['PDBchain'] = pdb_sifts_df['PDB']+pdb_sifts_df['CHAIN'] # add a column for PDBchain by joining

        Path('./Downloads').mkdir(parents=True, exist_ok=True)
        pdb_sifts_df.to_csv(Path('./Downloads/pdb_sifts.csv'), index=False)

def get_template_files():
    Template_list = []
    for template_res_num in [8, 7, 6, 5, 4, 3]:
        template_file = 'jess_templates_20230210/?_residues/results/?_res_templates.txt'.replace('?', str(template_res_num))

        with open(template_file, 'r') as f:
            lines = f.readlines()
        for line in lines:
            Template_list.append(line.strip())
            
    return Template_list

def evaluate_Temp_res_num():
    if not Path('./Downloads/Template_res_num_dict.json').is_file or not Path('./Downloads/Resnum_per_Template_dict.json').is_file():
        # Temp_res_num is not necessarily equal to the resnumber given in path:
        # Residues matching ANY amino acid type are not counted as they are too general
        # These have a value of 100 or hgiher in the second column indicating unspecific residues
        # Not all BB residues are unspecific! Some are targeted towards only Gly for example and thus have values < 100
        # Frankly these are still very promiscous and I'd like to exclude them too
        
        # Even if a Residue must match a particular AA, 6 atoms from a given residue may be selected
        # once for main chain and once for side chain
        # Therefore we only count unique template residues!
        
        # It seems that some templates even have the same 3-atom triplets at the same location twice. This I assume must be an error
        # again a reason to only count unique residues

        def count_resnums(Template):
            with open(Template, 'r') as f:
                Result = f.readlines()
            # We want to get a % of backbone RESIDUES in a template
            unique_residues = set()
            for line in Result:
                if line[:4] == 'ATOM':
                    if int(line[8:11]) < 100: # second column means it must be amino acid type specific
                        if line[17:20] != 'ANY': # exclude backbone atoms too - even if they are AA specific
                            res = int(line[22:26])
                            unique_residues.add(res)
            return len(unique_residues)

        Template_res_num_dict = {}
        Resnum_per_Template_dict = {}
        Resnum_per_Template_dict[3] = []
        Resnum_per_Template_dict[4] = []
        Resnum_per_Template_dict[5] = []
        Resnum_per_Template_dict[6] = []
        Resnum_per_Template_dict[7] = []
        Resnum_per_Template_dict[8] = []
        Template_list = get_template_files()
        with open(Path('./Downloads/excluded_templates.info'), 'w') as f:
            f.write('The following templates will not be evaluated due to insufficiently specified residues! \n')
            for Template in Template_list:
                residue_number = count_resnums(Template)
                Template_res_num_dict[Template] = residue_number
                if residue_number in range(3, 9):
                    Resnum_per_Template_dict[residue_number].append(Template)
                else:
                    f.write(Template + '\n')
        print('Downloads/excluded_templates.info')

        Path('./Downloads').mkdir(parents=True, exist_ok=True)
        with open(Path('./Downloads/Template_res_num_dict.json'), 'w') as f:
            json.dump(Template_res_num_dict, f)
        with open(Path('./Downloads/Resnum_per_Template_dict.json'), 'w') as f:
            json.dump(Resnum_per_Template_dict, f)
        
def template_relative_order():
    if not Path('./Downloads/Template_res_order_dict.json').is_file():
        # determine the relative sequence order of residues in a template

        Template_res_order_dict = {}
        Template_list = get_template_files()
        for Template in Template_list:
            res_ids = []
            with open(Template, 'r') as f:
                Result = f.readlines()
            for line in Result:
                if line[:4] == 'ATOM':
                    res_ids.append(int(line[22:26])) # get the residue ids to a list

            Residues = list(dict.fromkeys(res_ids)) # remove duplicates while preserving order
            res_rel_order = [sorted(Residues).index(i) for i in Residues] # determine the relative order of residues

            Template_res_order_dict[Template] = res_rel_order

        Path('./Downloads').mkdir(parents=True, exist_ok=True)
        with open(Path('./Downloads/Template_res_order_dict.json'), 'w') as f:
            json.dump(Template_res_order_dict, f)

def chunks(lst, n):
    """Yield successive n-sized chunks from lst."""
    for i in range(0, len(lst), n):
        yield lst[i:i + n]
            
def get_vec_template(residue):
    # dictionary in which the vectors from start to finish are defined for each residue type
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
    
    test_atom = list(residue)[0]
    res_type = test_atom[17:20].strip()
    vectup = vector_atom_type_dict[res_type]
    
    if vectup[1] != 'mid': # from first atom to second atom
        for index, atom in enumerate(residue):
            atom_name = atom[12:16].strip()
            if atom_name == vectup[0]:
                first_atom_index = index
                first_x = float(atom[30:38])
                first_y = float(atom[38:46])
                first_z = float(atom[46:54])
            elif atom_name == vectup[1]:
                second_atom_index = index
                second_x = float(atom[30:38])
                second_y = float(atom[38:46])
                second_z = float(atom[46:54])
                
        try:
            if first_x or second_x:
                pass
        except:
            return 'err', 'err', 'err'
                
        vec = [first_x-second_x, first_y-second_y, first_z-second_z]   
        return vec, [first_atom_index, second_atom_index], res_type
            
    else: # from the middle atom to the midpoint
        first_atom = False
        for index, atom in enumerate(residue):
            atom_name = atom[12:16].strip()
            if atom_name == vectup[0]:
                middle_atom_index = index
                middle_x = float(atom[30:38])
                middle_y = float(atom[38:46])
                middle_z = float(atom[46:54])
            elif first_atom == False:
                first_atom = True
                first_atom_index = index
                first_atom_name = atom_name
                side1_x = float(atom[30:38])
                side1_y = float(atom[38:46])
                side1_z = float(atom[46:54])
            elif first_atom == True:
                second_atom_index = index
                second_atom_name = atom_name
                side2_x = float(atom[30:38])
                side2_y = float(atom[38:46])
                side2_z = float(atom[46:54])
                
        try:
            if middle_x:
                pass
        except:
            return 'err', 'err', 'err'

        # calculate midpoint between the two side atoms        
        midpoint = [(side1_x+side2_x)/2, (side1_y+side2_y)/2, (side1_z+side2_z)/2]
        # calculate the vector between middle atom and midpoint
        vec = [middle_x-midpoint[0], middle_y-midpoint[1], middle_z-midpoint[2]]

        return vec, [middle_atom_index, first_atom_index, second_atom_index], res_type
        
def calculate_template_orientation_vectors(Template_list):
    Template_vec_dict = {}
    for Template in Template_list:
        if Template not in Template_vec_dict:
            Template_vec_dict[Template] = {}
        atom_lines = []
        with open(Template, 'r') as f: # read the Template
            Result = f.readlines()
            for line in Result:
                if line[:4] == 'ATOM':
                    atom_lines.append(line)

        residue_positions = [] # a list of positions for each atom in each residue (list of numpy arrays)
        residue_generator = chunks(atom_lines, 3)
        for index, residue in enumerate(residue_generator): # iterate over residues with index
            if index not in Template_vec_dict[Template]:
                Template_vec_dict[Template][index] = {}
            vec_all = get_vec_template(residue) # calculate the vector for each residue

            if vec_all[0] == 'err':
                print(Template)
                sys.exit()

            Template_vec_dict[Template][index]['vector'] = vec_all[0]
            Template_vec_dict[Template][index]['index_order'] = vec_all[1]
            Template_vec_dict[Template][index]['res_type'] = vec_all[2]

            res_pos = np.empty(shape=(3, 3)) # initialize an empty array
            for jndex, atom in enumerate(residue):
                atom_x = float(atom[30:38])
                atom_y = float(atom[38:46])
                atom_z = float(atom[46:54])
                res_pos[jndex] = np.array([atom_x, atom_y, atom_z]) # fill the array at a given row
            residue_positions.append(res_pos.tolist()) # to make it json serializable

            Template_vec_dict[Template]['residue_positions'] = residue_positions

    return Template_vec_dict
            
def check_other_files():
    Path('./data').mkdir(parents=True, exist_ok=True)

    EC_cofactor_mapping = Path('./data/EClist_cofactors_forRH.csv')
    AlphaFold_CATH = Path('./data/cath-v4_3_0.alphafold-v2.2022-11-22.tsv')
    PDBchain_to_CATH_Uniprot = Path('./data/pdb_chain_cath_uniprot.csv')
    PDBchain_to_EC_Uniprot = Path('./data/pdb_chain_enzyme.csv')
    CATH_names_mapping = Path('./data/cath-domain-list.txt')
    
    if not Path(EC_cofactor_mapping).is_file():
        sys.exit('EC to cofactor mapping csv from Neera is missing!')

    if not Path(PDBchain_to_CATH_Uniprot).is_file():
        urllib.request.urlretrieve('ftp://ftp.ebi.ac.uk/pub/databases/msd/sifts/flatfiles/csv/pdb_chain_cath_uniprot.csv.gz', './data/pdb_chain_cath_uniprot.csv.gz')
        
        with gzip.open('./data/pdb_chain_cath_uniprot.csv.gz', 'rb') as f_gz:
            with open('./data/pdb_chain_cath_uniprot.csv', 'wb') as f:
                shutil.copyfileobj(f_gz, f)
        # Mapping PDBchains to CATH name and Uniprot! Obtain via:
        # https://www.ebi.ac.uk/pdbe/docs/sifts/quick.html
                
        Path.unlink('./data/pdb_chain_cath_uniprot.csv.gz')

    if not Path(PDBchain_to_EC_Uniprot).is_file():

        urllib.request.urlretrieve('ftp://ftp.ebi.ac.uk/pub/databases/msd/sifts/flatfiles/csv/pdb_chain_enzyme.csv.gz', './data/pdb_chain_enzyme.csv.gz')
        
        with gzip.open('./data/pdb_chain_enzyme.csv.gz', 'rb') as f_gz:
            with open('./data/pdb_chain_enzyme.csv', 'wb') as f:
                shutil.copyfileobj(f_gz, f)

        Path.unlink('./data/pdb_chain_enzyme.csv.gz')

        # Mapping PDBchains to EC and Uniprot
        # https://www.ebi.ac.uk/pdbe/docs/sifts/quick.html

    if not Path(CATH_names_mapping).is_file():
        # Mapping CATH ids to CATH numbers and CATH domain names
        urllib.request.urlretrieve('ftp://orengoftp.biochem.ucl.ac.uk/cath/releases/latest-release/cath-classification-data/cath-domain-list.txt', './data/cath-domain-list.txt')
    if not Path(AlphaFold_CATH).is_file():
        # Mapping UniProt AFDB Proteins to CATH Domains
        urllib.request.urlretrieve('ftp://orengoftp.biochem.ucl.ac.uk/alphafold/cath-v4.3.0-model-organisms/cath-v4_3_0.alphafold-v2.2022-11-22.tsv', './data/cath-v4_3_0.alphafold-v2.2022-11-22.tsv')
    
def main():
    get_mcsa_ec()
    get_mcsa_cath()
    get_mcsa_interpro()
    evaluate_Temp_res_num()
    template_relative_order()
    calculate_template_orientation_vectors()
    check_other_files()
    make_pdb_sifts_df()
    
    print('All files for annotation of jess results are there!')
        
if __name__ == "__main__":
    main()
    
