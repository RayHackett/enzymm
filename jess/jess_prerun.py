"""
written by Ray
jess_prerun.py

Generate and or check files necessary to run jess:

Outsources to:
    common_functions.py

Inputs:
    files with paths to each each template split by template residue number
    
Outputs:
    /homes/hackett/Downloads/MCSA_EC_mapping.json
    /homes/hackett/Downloads/MCSA_CATH_mapping.json
    /homes/hackett/Downloads/MCSA_interpro_dict.json
    /homes/hackett/Downloads/Template_res_num_dict.json
    /homes/hackett/Downloads/Resnum_per_Template_dict.json
    
"""

import sys
from pathlib import Path
import json
import re
import numpy as np
import pandas as pd
from Bio.PDB import PDBParser
from Bio.PDB import MMCIFParser
from Bio.PDB.DSSP import DSSP

"""
Comments on Biopython - Version modified by Yannis:
# changes to DSSP.py see https://swift.cmbi.umcn.nl/gv/dssp/DSSP_3.html
# if chainid (line[11] == '>' then look at the end of the line)

#MMCIF parser assumes that models are arranged contionously - they are not example: 1iv4
#it should instead return to the actually referenced model in the structure instead of incrementing by 1
# we bypass this by adding in lines:

491                 try: 
492                     chain = model[chain_id]
493                 except:
494                     print('DSSP warning in file', in_file)
495                     break

Note that this is not a fix!!!
"""

d = Path(__file__).resolve().parent

from ..common_functions import request_url

def get_mcsa_ec():
    # m-csa entry map to EC number
    # every m-csa entry only has one EC number
    if not Path(d, '../Downloads/MCSA_EC_mapping.json').is_file():
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
                    
            with open(Path(d, '../Downloads/MCSA_EC_mapping.json'), 'w') as f:
                json.dump(ec_dict, f)

def get_mcsa_cath():
    # m-csa entry map to CATH id
    # entry have have multiple CATH ids
    if not Path(d, '../MCSA_CATH_mapping.json').is_file():
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
                
            with open(Path(d, '../MCSA_CATH_mapping.json'), 'w') as f:
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
    if not Path(d, '../MCSA_interpro_dict.json').is_file():
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
        
        with open(Path(d, '../Downloads/MCSA_interpro_dict.json'), 'w') as f:
                json.dump(MCSA_interpro_dict, f)
                
def make_pdb_sifts_df():
    if not Path(d, '../Downloads/pdb_sifts.csv').is_file():
        PDBchain_to_CATH_Uniprot = Path(d, '../Downloads/pdb_chain_cath_uniprot.csv')
        PDBchain_to_EC_Uniprot = Path(d, '../Downloads/pdb_chain_enzyme.csv')
        CATH_names_mapping = Path(d, '../Downloads/cath-domain-list.txt')

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
        pdb_sifts_df.to_csv(Path(d, '../Downloads/pdb_sifts.csv'), index=False)

def get_template_files():
    Template_list = []
    for template_res_num in [8, 7, 6, 5, 4, 3]:
        template_file = '/hps/nobackup/thornton/hackett/jess_templates_20230210/?_residues/results/?_res_templates.txt'.replace('?', str(template_res_num))

        with open(template_file, 'r') as f:
            lines = f.readlines()
        for line in lines:
            Template_list.append(line.strip())
            
    return Template_list

def evaluate_Temp_res_num():
    if not Path(d, '../Downloads/Template_res_num_dict.json').is_file or not Path(d, '../Downloads/Resnum_per_Template_dict.json').is_file():
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
        with open(Path(d, '../Downloads/excluded_templates.info'), 'w') as f:
            f.write('The following templates will not be evaluated due to insufficiently specified residues! \n')
            for Template in Template_list:
                residue_number = count_resnums(Template)
                Template_res_num_dict[Template] = residue_number
                if residue_number in range(3, 9):
                    Resnum_per_Template_dict[residue_number].append(Template)
                else:
                    f.write(Template + '\n')
        print('See ../Downloads/excluded_templates.info')

        with open(Path(d, '../Downloads/Template_res_num_dict.json'), 'w') as f:
            json.dump(Template_res_num_dict, f)
        with open(Path(d, '../Downloads/Resnum_per_Template_dict.json'), 'w') as f:
            json.dump(Resnum_per_Template_dict, f)
        
def template_relative_order():
    if not Path(d, '../Downloads/Template_res_order_dict.json').is_file():
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

        with open(Path(d, '../Downloads/Template_res_order_dict.json'), 'w') as f:
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
        
def calculate_template_orientation_vectors():
    if not Path(d, '../Downloads/Template_vec_dict.json').is_file():
        Template_vec_dict = {}
        Template_list = get_template_files() # load all the Template files into a list
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

        with open(Path(d, '../Downloads/Template_vec_dict.json'), 'w') as f:
            json.dump(Template_vec_dict, f)
            
def calculate_template_dssp():
    # Calculates Solvent Accessible Surface Area for all Templates by residue and stores the info in dictionaries
    def get_dssp(pdb):
        # see https://swift.cmbi.umcn.nl/gv/dssp/DSSP_3.html
        dssp_path = '/hps/software/users/thornton/hackett/anaconda3/envs/mod_biopy/bin/mkdssp'
        auth_structure = MMCIFParser(auth_chains=True, auth_residues=True, QUIET=True).get_structure(str(Path(pdb).name), Path(pdb))
        #label_structure = MMCIFParser(auth_chains=False, auth_residues=True, QUIET=True).get_structure(str(Path(pdb).name), Path(pdb))
        model = auth_structure[0] # always only one model per pdb file since we are dealing with experimental Non-NMR or predicted structures
        try:
            dssp = DSSP(model, Path(pdb), dssp=dssp_path, acc_array='Miller') # Miller since the Group useses that
        except:
            dssp = False
        return auth_structure, dssp  
    
    if not Path(d, '../Downloads/Template_dssp_dict.json').is_file():
        # Note that I did not parallelise this as I really only have to run it once...
        print('Calculating DSSP for all Template Structures. This might take a while ...')
        
        Template_list = get_template_files() # load all the Template files into a list
        # get all the pdb identifiers
        template_pdbs = set()
        for Template in Template_list:
            template_pdbs.add(Path(Template.strip()).name.split('.')[2][:4])
        # get a list of the cif assemblies for each pdb identifier
        pdb_assembly_mapping = []
        with open('/nfs/research/thornton/riziotis/research/phd/data/pdbe/assembly_cif/cif.ls', 'r') as f:
            cif_files = f.readlines()
        for cif_file in cif_files:
            if cif_file[4:8] in template_pdbs:
                pdb_assembly_mapping.append('/nfs/research/thornton/riziotis/research/phd/data/pdbe/assembly_cif/' + cif_file[:19] + '.cif')
                
        # calculate DSSP for all the Template CIF files and store the info in a dictionary
        template_dssp_dict = {}
        failed_dssp_list = []
        for pdb in pdb_assembly_mapping:
            if Path(pdb):
                auth_structure, dssp = get_dssp(pdb)
                if not dssp: # if Dssp could not be calculated, skip it and print the name of the structure
                    failed_dssp_list.append(pdb)
                    continue
                    
                # convert to dictionary. this fails for multi-model NMR structures
                try:
                    dssp_dict = dict(dssp) # keys are (chainid, (' ', Exp_ResNumber, ' '))
                except:
                    failed_dssp_list.append(pdb)
                    continue
                # value is a tuple of:
                # (dssp index, amino acid, secondary structure, relative ASA, phi, psi,
                # NH_O_1_relidx, NH_O_1_energy, O_NH_1_relidx, O_NH_1_energy,
                # NH_O_2_relidx, NH_O_2_energy, O_NH_2_relidx, O_NH_2_energy)

                pdb_name = str(Path(pdb).name.split('.')[0][:4])
                if pdb_name not in template_dssp_dict:
                    template_dssp_dict[pdb_name] = {}

                rsa_list = []
                for i, residue in enumerate(auth_structure.get_residues()):
                    # check if it is a value
                    data = dssp_dict.get((residue.get_full_id()[2], (' ', residue.id[1], residue.id[2])))
                    if data and isinstance(data[3],float): # incomplete residues are not evaluated, HETATM get 'NA'
                        rsa = float(data[3])
                        #ss = data[2] # secondary structure
                    else:
                        rsa = 9.99
                    rsa_list.append(rsa)
                    #plddt = residue['CA'].get_bfactor() / 100.0 # plddt

                # yannis uses the auth residue numbers in his templates
                for i, residue in enumerate(auth_structure.get_residues()):
                    resname = str(residue.get_full_id()[2]+str(residue.id[1]))
                    template_dssp_dict[pdb_name][resname] = (rsa_list[i], residue.get_resname())
        
        # DSSP will fail on the following NMR structures! 
        # Dictionary with file name as key and list of catalytic residues as values (with chain ID)
        known_failures = {'1i4v-assembly-1.cif': ['A60','A93','A59','A97'],
                          '1qfn-assembly-1.cif': ['A13','A11','A18','A72','A8','A10','A14'],
                          '2xy8-assembly-1.cif': ['A167','A162','A61','A12','A14']}
        
        def add_manual_dssp(file, catalytic_residues):
            template_dssp_dict[Path(file).name[:4]] = {}
            for i in catalytic_residues:
                auth_structure, dssp = get_dssp(file)
                resname = i
                number = int(i[1:])
                aa_code = auth_structure[0]["A"][number].get_resname()
                template_dssp_dict[Path(file).name[:4]][resname] = (float(dssp[('A', (' ', number, ' '))][3]), aa_code)
        
        # add the kown failures manually. if one is missing, add it to the known_failures dictionary
        for failure in failed_dssp_list:
            name = Path(failure).name 
            if name in known_failures.keys():
                add_manual_dssp(failure, known_failures[name])
            else:
                sys.exit("""DSSP failed on a structure for which no catalytic residues were given.\n
                The structure was {}.\n
                Please fix this""".format(failure))
                
        with open(Path(d, '../Downloads/Template_dssp_dict.json'), 'w') as f:
            json.dump(template_dssp_dict, f)
            
    # now we want to cut down on space needed by only saving the DSSP values for the residues which are actually in the template files themselves
    # make a dict with template_name, residue_name and then the solvation value
    if not Path(d, '../Downloads/Template_res_dssp_dict.json').is_file():
        # This dictionary only includes those residues which are in the template file itself - so only the catalytic ones
        template_res_dssp_dict = {}
        baddies = set()
        Template_list = get_template_files()
        for Template in Template_list:
            with open(Template, 'r') as f:
                template_lines = f.readlines()
            pdb_name = Path(Template).name.split('.')[2][:4]

            template_res_dssp_dict[Template] = {}

            atom_lines = []
            for line in template_lines:
                if line[:4] == 'ATOM':
                    atom_lines.append(line)

            for line in atom_lines[::3]:
                chain = line[20:22].strip()
                resname = line[17:20]
                res_id = int(line[22:26])
                temp_atom = chain + str(res_id)
                # only add the solvation value if the residue name matches. this is a small check. it is a proof of correct mapping
                try:
                    #template_res_dssp_dict[Template][temp_atom] = new_dict[pdb_name][temp_atom][0]
                    if resname == 'ANY':
                        template_res_dssp_dict[Template][temp_atom] = template_dssp_dict[pdb_name][temp_atom][0]
                    elif resname == template_dssp_dict[pdb_name][temp_atom][1]: # cannonical amino acids
                        template_res_dssp_dict[Template][temp_atom] = template_dssp_dict[pdb_name][temp_atom][0]
                    elif resname == 'PTM': # PTM gets added if the reference is a PTM - not necessarily the structure from which the template was derived
                        template_res_dssp_dict[Template][temp_atom] = template_dssp_dict[pdb_name][temp_atom][0]
                    else:
                        print('resnames dont match')
                        print(Template, temp_atom)
                except:
                    baddies.add(pdb_name)
        if baddies:
            sys.exit("""DSSP annotations failed for the following structures as residue names were not mapped correctly: \n
            {}""".format(baddies))

        with open(Path(d, '../Downloads/Template_res_dssp_dict.json'), 'w') as f:
            json.dump(template_res_dssp_dict, f)
            
def check_other_files():
    EC_cofactor_mapping = Path(d, '../Downloads/EClist_cofactors_forRH.csv')
    AlphaFold_CATH = Path(d, '../Downloads/cath-v4_3_0.alphafold-v2.2022-11-22.tsv')
    PDBchain_to_CATH_Uniprot = Path(d, '../Downloads/pdb_chain_cath_uniprot.csv')
    PDBchain_to_EC_Uniprot = Path(d, '../Downloads/pdb_chain_enzyme.csv')
    CATH_names_mapping = Path(d, '../Downloads/cath-domain-list.txt')
    
    if not Path(EC_cofactor_mapping).is_file():
        sys.exit('EC to cofactor mapping csv from Neera is missing!')
    if not Path(PDBchain_to_CATH_Uniprot).is_file():
        sys.exit("""Mapping PDBchains to CATH name and Uniprot is missing! Obtain via: \n
        https://www.ebi.ac.uk/pdbe/docs/sifts/quick.html \n
        ftp://ftp.ebi.ac.uk/pub/databases/msd/sifts/flatfiles/csv/pdb_chain_cath_uniprot.csv.gz""")
    if not Path(PDBchain_to_EC_Uniprot).is_file():
        sys.exit(""" Mapping PDBchains to EC and Uniprot is missing! Obtain via: \n
        https://www.ebi.ac.uk/pdbe/docs/sifts/quick.html \n
        ftp://ftp.ebi.ac.uk/pub/databases/msd/sifts/flatfiles/csv/pdb_chain_enzyme.csv.gz""")
    if not Path(CATH_names_mapping).is_file():
        sys.exit("""Mapping CATH ids to CATH numbers and CATH domain names is missing! \n
        Obtain through ftp://orengoftp.biochem.ucl.ac.uk/cath/releases/latest-release/cath-classification-data/cath-domain-list.txt""")
    if not Path(AlphaFold_CATH).is_file():
        sys.exit("""Alphafold CATH mapping is missing. Obtain via: \n
            ftp ftp://orengoftp.biochem.ucl.ac.uk/
            /alphafold/cath-v4.3.0-model-organisms/cath-v4_3_0.alphafold-v2.2022-11-22.tsv""")
    
def main():
    get_mcsa_ec()
    get_mcsa_cath()
    get_mcsa_interpro()
    evaluate_Temp_res_num()
    template_relative_order()
    calculate_template_orientation_vectors()
    calculate_template_dssp()
    check_other_files()
    make_pdb_sifts_df()
    
    print('All files for annotation of jess results are there!')
        
if __name__ == "__main__":
    main()
    
