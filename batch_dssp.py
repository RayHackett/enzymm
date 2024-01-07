"""
written by Ray
script 1: batch_dssp.py
script to split an input list of pdb or cif files into managable batches which are then submitted to the dssp script
outputs on if the dssp script ran successfully for each structure are then collected and aggregated

Outsources to:
    common_functions.py
    LSF (bsub)
    dssp_run.py

Inputs:
    list of pdb structure files supplied with -i
    batch sizes supplied with -s
    
Outputs:
    ./collected_dssp_successes.list with all the successfully dssp scored structures

"""

import sys
import re
import time
from pathlib import Path
import subprocess
import argparse
import uuid

from common_functions import file_batching

parser = argparse.ArgumentParser()
parser.add_argument('-i', '--input', help='file with list of protein structures to scan')
parser.add_argument('-s', '--size', help='Batch size for how many structures to run jess against')
args = parser.parse_args()

pdb_file = args.input
batch_size = args.size

cwd = Path.cwd()

results_path = Path(cwd,'dssp_results') 
results_path.mkdir(exist_ok=True)

########### Split pdbs into Batches ################

# pdb_file comes from argparse -i, batch_size comes from -s
batch_list = file_batching(pdb_file, int(batch_size), results_path, 'dssp_batch')

########### submit jobs ############################
unique_tag = str(uuid.uuid4())

for file_name in batch_list:
    batch_number = re.findall(r'\d+', str(file_name))[-1]
    job_name = 'dssp_batch_' + unique_tag + batch_number
    err_name = Path(results_path,'dssp_log' + str(batch_number) + '.err')
    out_name = Path(results_path,'dssp_log' + str(batch_number) + '.out')
    dssp_run_script = './dssp_run.py'
    script = ' '.join(['bsub', '-J', job_name, '-e', str(err_name), '-o', str(out_name), 
                            "'", 'python', dssp_run_script, '-i', str(file_name), "'"])   
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
        with open('failed_dssp_job_submissions.txt', 'a') as f:
            f.write(file_name + '\n')
                
##########  Wait for all jobs to finish  ###########
try:
    subprocess.check_call("bwait -w 'ended({})'".format('dssp_batch_' + unique_tag + '*'), shell=True)
except:
    pass

############# Concatenate Output #####################
file_template = 'dssp_batch?.list' 
batches = list(range(1,len(batch_list)+1))

all_lines = []
for batch_number in batches:
    file = file_template.replace('?', str(batch_number))
    file_name = Path(results_path, file)
    if file_name.is_file():
        with open(file_name) as f:
            all_lines.extend(f.readlines())

with open('collected_dssp_successes.list', 'w') as f:
    for i in all_lines:
        f.write(i)
            
