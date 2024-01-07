"""
written by Ray
script: jess_submission.py
splits pdb_file into batches of a given size and submits each batch to the jess pipeline in jess_run.py

Outsources to:
    common_functions.py
    LSF (bsub)
    jess_prerun.py
    jess_run.py

Inputs:
    Reads list of pdb files supplied after -i flag
    pdb Batch size for jess calculations after -s flag
    Jess parameters rmsd, distance, max_dynamic_distance score_cutoff after -j flag
    
Outputs:
    ./jess/pdb_batch<number>.txt file with batch of pdb files
    ./jess results from jess seperated by residue number of the templates
    ./jess_logs log files from jess runs
    ./<batch>all_res_info.json file with all info and annotations for all hits to the queries in the batch
    ./combined_all_res_info.json
    ./jess.done
"""

import sys
import subprocess
import re
import time
import json
import argparse
from pathlib import Path
import uuid

sys.path.append('/hps/software/users/thornton/hackett/software')
from common_functions import file_batching
sys.path.append('/hps/software/users/thornton/hackett/jess')
import jess_prerun

def main(pdb_file, batch_size, jess_params):
    
    ############# Setting paths and making directories ###########
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
    
    # jess log path
    global jess_log_path
    jess_log_path = Path(cwd,'jess_logs') 
    jess_log_path.mkdir(exist_ok=True)
    
    ####### jess_prerun.py
    jess_prerun.get_mcsa_ec()
    jess_prerun.get_mcsa_cath()
    jess_prerun.get_mcsa_interpro() # takes a bit longer too!
    jess_prerun.evaluate_Temp_res_num()
    jess_prerun.template_relative_order()
    jess_prerun.calculate_template_orientation_vectors() # a bit expensive!
    jess_prerun.calculate_template_dssp() # takes about 20min!
    jess_prerun.check_other_files()
    jess_prerun.make_pdb_sifts_df()
    
    ########### Split pdbs into Batches ################
    
    # pdb_file comes from argparse -i, batch_size comes from -s
    batch_list = file_batching(pdb_file, int(batch_size), jess_path, 'pdb_batch')
    
    ######## wirte the template files ###############
    
    # load dictionary mapping residue_number to Templates
    # These are the Templates which are then read by Jess
    with open('/homes/hackett/Downloads/Resnum_per_Template_dict.json', 'r') as f:
        Resnum_per_Template_dict = json.load(f)
    # json always converts keys to strings. we want int type
    Resnum_per_Template_dict = {int(k): v for k,v in Resnum_per_Template_dict.items()}
    
    for template_res_num in [8, 7, 6, 5, 4, 3]:
        template_file = Path(template_folder, 'list_of_{}_residue_templates.list'.format(str(template_res_num)))
        with open(template_file, 'w') as f:
            list_of_templates = Resnum_per_Template_dict[template_res_num]
            for template in list_of_templates:
                # get only templates of cluster 1
                # get the first character after cluster_
                #if template.split('cluster_')[1][0] == '1':
                    #f.write(template + '\n')
                f.write(template + '\n')
                
    ######### Submit Jess run jobs #################
    
    unique_tag = str(uuid.uuid4())
    
    for file_name in batch_list:
        batch_number = re.findall(r'\d+', str(file_name))[-1]
        job_name = 'jess_batch_' + unique_tag + batch_number
        err_name = Path(jess_log_path,'jess_log' + str(batch_number) + '.err')
        out_name = Path(jess_log_path,'jess_log' + str(batch_number) + '.out')
        jess_run_script = '/hps/software/users/thornton/hackett/jess/jess_run.py'
        script = ' '.join(['bsub', '-J', job_name, '-e', str(err_name), '-o', str(out_name), 
                            "'", 'python', jess_run_script, '-i', str(file_name), '-j', jess_params, "'"])
        
        attempts = 0
        job_status = False
        while job_status == False and attempts < 3:
            try:
                subprocess.check_call(script, shell=True)
                job_status = True
            except subprocess.CalledProcessError as e:
                attempts += 1
                time.sleep(3)
                print(e)
                print('retrying job submission in 3s...')
        if attempts == 3:
            with open('failed_jess_job_submissions.txt', 'w') as f:
                f.write(file_name + '\n')
                
    ##########  Wait for all jobs to finish  ###########
    try:
        subprocess.check_call("bwait -w 'ended({})'".format('jess_batch_' + unique_tag + '*'), shell=True)
    except:
        pass
    
    ############# Concatenate Output #####################
    # might have to change this pattern if the working_file in jess_run changes
    file_template = 'pdb_batch?all_res_info.json' 
    batches = list(range(1,len(batch_list)+1))

    all_res_dict = {}
    for batch_number in batches:
        file = file_template.replace('?', str(batch_number))
        file_name = Path(cwd,'jess_results/' + file)
        if file_name.is_file():
            with open(file_name) as f:
                new_dict = json.load(f)

            for key, value_keys in new_dict.items():
                if key not in all_res_dict:
                    all_res_dict[key] = {}
                for sub_key, match_list in value_keys.items():
                    if sub_key not in all_res_dict[key]:
                        all_res_dict[key][sub_key] = []
                    for match in match_list:
                        all_res_dict[key][sub_key].append(match)
                
    with open('combined_all_res_info.json', 'w') as f:
            json.dump(all_res_dict, f)
    
    ################### Finished  #################################
    message = 'Jess Calculations finished succesfully'
    
    with open ('jess.done', 'w') as f:
        f.write(message)
        
if __name__ == "__main__":
    
    parser = argparse.ArgumentParser()
    parser.add_argument('-i', '--input', help='file with list of protein structures to scan')
    parser.add_argument('-s', '--size', help='Batch size for how many structures to run jess against')
    parser.add_argument('-j', '--jess', nargs = '+', help='Jess space seperated parameters rmsd, distance, max_dynamic_distance, score_cutoff, optional flags as a string')
    args = parser.parse_args()
    
    pdb_file = args.input
    batch_size = args.size
    jess_params = ' '.join(args.jess)
    
    main(pdb_file, batch_size, jess_params)

