#!/bin/sh
#SBATCH -t 24:00:00
#SBATCH -e Swissprot_%a_%j.err
#SBATCH -o Swissprot_%a_%j.out
#SBATCH --qos=normal
#SBATCH --mem=128G
#SBATCH --nodes=1
#SBATCH --no-requeue

# set bash strict mode http://redsymbol.net/articles/unofficial-bash-strict-mode/
set -euo pipefail
IFS=$'\n\t'

TIMESTAMP=$(date "+%Y%m%d_%H%M%S")
PROFILE="cluster_lumc" # local or cluster_lumc or cluster_embl
PROJECT="template_matching/Human_Swissprot"

if [ ${PROFILE} == "cluster_lumc" ]; then
    WORKFLOW="/exports/archive/lucid-grpzeller-primary/hackett/template_matching/nextflow/template_matcher.nf"
    CONFIG="/exports/archive/lucid-grpzeller-primary/hackett/template_matching/nextflow/matcher.config"
    module load bioinformatics/tools/Nextflow/22.10.6 # LUMC cluster
    WORKDIR=/exports/lucid-grpzeller-work/rehackett/${PROJECT}/
    PARAMS="params_lumc.yml"

elif [ ${PROFILE} == "cluster_embl" ]; then
    module load Nextflow/23.10.1 # EMBL cluster
    WORKDIR=/scratch/rhackett/${PROJECT}/
    PARAMS="params_embl.yml"

elif [ ${PROFILE} == "local" ]; then
    mamba activate matcher
    WORKDIR=~/Documents/${PROJECT}/
    PARAMS="params.yml"
fi

## will be created
mkdir -p ${WORKDIR}
## create a symlink from current working dir to this location
# ln -sf $WORKDIR work

# nextflow pull ${WORKFLOW} -r ${WORKFLOW_VERSION}

nextflow run ${WORKFLOW} \
    -profile ${PROFILE} \
    -c ${CONFIG} \
    -params-file ${PARAMS} \
    -work-dir ${WORKDIR} \
    -with-conda \
    -with-trace trace.${TIMESTAMP}.txt \
    -with-report report.${TIMESTAMP}.html \
    -resume

# 	-r ${WORKFLOW_VERSION} \ # add this if working with remote git pipeline versions
