#!/bin/bash

# Script to download files from RCSB http file download services.
# Use the -h switch to get help on usage.
# See https://www.rcsb.org/docs/programmatic-access/batch-downloads-with-shell-script
# modified to download specific assemblies!
# specify assembly like so: {pdb_id}-{assembly_number}
# assumes a 1 digit assembly number btw

if ! command -v curl &> /dev/null
then
    echo "'curl' could not be found. You need to install 'curl' for this script to work."
    exit 1
fi

PROGNAME=$0
BASE_URL="https://files.rcsb.org/download"

usage() {
    cat << EOF >&2
Usage: $PROGNAME -f <file> [-o <dir>] [-c] [-p]

    -f <file>: the input file containing a comma-separated list of PDB ids
    -o  <dir>: the output dir, default: current dir
    -c       : download a cif.gz file for each PDB id
    -p       : download a pdb.gz file for each PDB id (not available for large structures)
    -a       : download a pdbX.gz file (Xth bioassembly) for each PDB id (not available for large structures)
    -A       : download an assemblyX.cif.gz file (Xth bioassembly) for each PDB id
    -x       : download a xml.gz file for each PDB id
    -s       : download a sf.cif.gz file for each PDB id (diffraction only)
    -m       : download a mr.gz file for each PDB id (NMR only)
    -r       : download a mr.str.gz for each PDB id (NMR only)
EOF
    exit 1
}

download() {
    url="$BASE_URL/$1"
    out="$2/$1"
    
    [ -f "$out" ] && return
    
    echo "Downloading $url to $out"
    curl -s -f "$url" -o "$out" || echo "Failed to download $url"                 
}  

listfile=""
outdir="."
cif=false
pdb=false
pdb1=false
cifassembly1=false
xml=false
sf=false
mr=false
mrstr=false
while getopts f:o:cpaAxsmr o
do
    case $o in
        (f) listfile=$OPTARG;;
        (o) outdir=$OPTARG;;
        (c) cif=true;;
        (p) pdb=true;;
        (a) pdb1=true;;
        (A) cifassembly1=true;;
        (x) xml=true;;
        (s) sf=true;;
        (m) mr=true;;
        (r) mrstr=true;;
        (*) usage
    esac
done
shift "$((OPTIND - 1))"

if [ "$listfile" == "" ]
then
    echo "Parameter -f must be provided"
    exit 1
fi
contents=$(cat $listfile)

# see https://stackoverflow.com/questions/918886/how-do-i-split-a-string-on-a-delimiter-in-bash#tab-top
IFS=',' read -ra tokens <<< "$contents"

for token in "${tokens[@]}"
do
    IFS='-' read -ra assembly_id <<< "$token"

    if [ "$cif" == true ]
    then
        download ${assembly_id[0]}.cif.gz $outdir
    fi
    if [ "$pdb" == true ]
    then
        download ${assembly_id[0]}.pdb.gz $outdir
    fi
    if [ "$pdb1" == true ]
    then
        download ${assembly_id[0]}.pdb${assembly_id[1]}.gz $outdir
    fi
    if [ "$cifassembly1" == true ]
    then
        download ${assembly_id[0]}-assembly${assembly_id[1]}.cif.gz $outdir
    fi
    if [ "$xml" == true ]
    then
        download ${assembly_id[0]}.xml.gz $outdir
    fi
    if [ "$sf" == true ]
    then
        download ${assembly_id[0]}-sf.cif.gz $outdir
    fi
    if [ "$mr" == true ]
    then
        download ${assembly_id[0]}.mr.gz $outdir
    fi
    if [ "$mrstr" == true ]
    then
        download ${assembly_id[0]}_mr.str.gz $outdir
    fi

done







