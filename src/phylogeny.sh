#!/bin/bash

if [ $# -lt 2 ]
then
    echo ""
    echo "Usage: $0 <outprefix> <sample1.cons.fa> <sample2.cons.fa> ... <sampleN.cons.fa>"
    echo ""
    exit -1
fi

SCRIPT=$(readlink -f "$0")
BASEDIR=$(dirname "$SCRIPT")
export PATH=${BASEDIR}/../conda/bin:${PATH}

# Input parameters
OUTP=${1}
shift 1
REF=${BASEDIR}/../ref/NC_045512.2.fa
THREADS=4

# Compute MSA
cat $@ | gzip -c > ${OUTP}.fa.gz
mafft --auto --thread ${THREADS} --keeplength --addfragments <(zcat ${OUTP}.fa.gz) ${REF} > ${OUTP}.msa

# Phylogenetic inference
iqtree -nt ${THREADS} -s ${OUTP}.msa

# Plotting
Rscript ${BASEDIR}/../scripts/phylogeny.R ${OUTP}.msa.treefile 
