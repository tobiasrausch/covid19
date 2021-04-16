#!/bin/bash

if [ $# -ne 4 ]
then
    echo ""
    echo "Usage: $0 <global|europe> <outprefix> <multi.meta.tsv> <multi.cons.fa>"
    echo ""
    exit -1
fi

SCRIPT=$(readlink -f "$0")
BASEDIR=$(dirname "$SCRIPT")
export PATH=${BASEDIR}/../conda/bin:${PATH}

# Input parameters
REGION=${1}
OUTP=${2}
DATA=${BASEDIR}/../../../covid/nextstrain/
THREADS=4

# Concatenate data
cp ${DATA}/ncov_${REGION}.fasta ${BASEDIR}/../ncov/data/example_sequences.fasta
cp ${DATA}/ncov_${REGION}.tsv ${BASEDIR}/../ncov/data/example_metadata.tsv
cat ${3} >> ${BASEDIR}/../ncov/data/example_metadata.tsv
cat ${4} >> ${BASEDIR}/../ncov/data/example_sequences.fasta

# Run ncov
cd ${BASEDIR}/../ncov/
rm -rf auspice benchmarks logs results
export AUGUR_RECURSION_LIMIT=10000
snakemake --cores 4 --profile ./my_profiles/getting_started

