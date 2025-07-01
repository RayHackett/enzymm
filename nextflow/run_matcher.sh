#!/bin/bash

#SBATCH -t 01:00:00
#SBATCH -e NM_%a_%j_err.txt
#SBATCH -o NM_%a_$j_out.txt
#SBATCH --qos=normal
#SBATCH --mem 8G
#SBATCH --no-requeue

# set bash strict mode http://redsymbol.net/articles/unofficial-bash-strict-mode/
set -euo pipefail
IFS=$'\n\t'

TIMESTAMP=$(date "+%Y%m%d_%H%M%S")
PROFILE="cluster" # local or cluster
PROJECT="testing_nextflow"

# # use with github
# WORKFLOW=cschu/nevermore
# WORKFLOW_VERSION=main

## my local version
WORKFLOW="/exports/archive/lucid-grpzeller-primary/hackett/EBI_Project/src/nextflow/template_matcher.nf"

if [ ${PROFILE} == "cluster" ]; then
    module load bioinformatics/tools/Nextflow/22.10.6 # LUMC cluster

elif [ ${PROFILE} == "local" ]; then
    mamba activate matcher

fi

## will be created
WORKDIR=/exports/lucid-grpzeller-work/hackett/${PROJECT}/
mkdir -p ${WORKDIR}
## create a symlink from current working dir to this location
ln -sf $WORKDIR work

# nextflow pull ${WORKFLOW} -r ${WORKFLOW_VERSION}

nextflow run ${WORKFLOW} \
    -profile ${PROFILE} \
	-c matcher.config \
	-params-file params.yml \
	-work-dir ${WORKDIR} \
	-with-trace trace.${TIMESTAMP}.txt \
	-with-report report.${TIMESTAMP}.html \
	-resume

# 	-r ${WORKFLOW_VERSION} \ # add this if working with remote git pipeline versions
