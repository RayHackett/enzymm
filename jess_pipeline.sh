####################################################################################################################
# script 1: batch_dssp.py
# script to split an input list of pdb or cif files into managable batches which are then submitted to the dssp script
# outputs on if the dssp script ran successfully for each structure are then collected and aggregated
#
# Outsources to:
#    common_functions.py
#     LSF (bsub)
#     dssp_run.py
#
# Inputs:
#     list of pdb structure files supplied with -i 
#     batch sizes supplied with -s
#     
# Outputs:
#     ./collected_dssp_successes.list with all the successfully dssp scored structures
# 

if [ -f "collected_dssp_successes.list" ]; then
    if [ "collected_dssp_successes.list" -ot "./succesful_cifs.list" ]; then
        python /hps/software/users/thornton/hackett/pipeline/batch_dssp.py -i succesful_cifs.list -s 100
    fi
else
    python /hps/software/users/thornton/hackett/pipeline/batch_dssp.py -i succesful_cifs.list -s 100
fi

####################################################################################################################
# written by Ray
# script: jess_submission.py
# splits pdb_file into batches of a given size and submits each batch to the jess pipeline in jess_run.py
# 
# Outsources to:
#     common_functions.py
#     LSF (bsub)
#     jess_prerun.py
#     jess_run.py
# 
# Inputs:
#     Reads list of pdb files supplied after -i flag
#     pdb Batch size for jess calculations after -s flag
#     Jess parameters rmsd, distance, max_dynamic_distance score_cutoff after -j flag
#     optional flags as a string like n: do not transform coordinates of hit into	the template coordinate frame
#     others are fiqe (see jess for description)
#     
# Outputs:
#     ./jess/pdb_batch<number>.txt file with batch of pdb files
#     ./jess results from jess seperated by residue number of the templates
#     ./jess_logs log files from jess runs
#     ./<batch>all_res_info.json file with all info and annotations for all hits to the queries in the batch
#     ./combined_all_res_info.json
#     ./jess.done

if [ -f "jess.done" ]; then
    if [ "jess.done" -ot "./pdbs/collected_dssp_successes.list" ]; then
        python /hps/software/users/thornton/hackett/pipeline/jess_submission.py -i collected_dssp_successes.list -s 2 -j 2 1.0 1.0 0
    fi
else
    python /hps/software/users/thornton/hackett/pipeline/jess_submission.py -i collected_dssp_successes.list -s 2 -j 2 1.0 1.0 0
fi

################################################################################################################

