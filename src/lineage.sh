#!/bin/bash

if [ $# -ne 2 ]
then
    echo ""
    echo "Usage: $0 <input.fasta> <unique_sample_id>"
    echo ""
    exit -1
fi

SCRIPT=$(readlink -f "$0")
BASEDIR=$(dirname "$SCRIPT")
export PATH=${BASEDIR}/../conda/bin:${PATH}

# Activate pangolin environment
source activate pangolin

# Input parameters
FASTA=${1}
OUTP=${2}
THREADS=4

# Run pangolin
if [ -f ${FASTA} ]
then
    pangolin -t ${THREADS} --outfile ${OUTP}.lineage.csv ${FASTA}
fi

# Deactivate
conda deactivate
