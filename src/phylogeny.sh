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

source activate	covid19

# Input parameters
OUTP=${1}
shift 1
REF=${BASEDIR}/../ref/NC_045512.2.fa
THREADS=4

# Concatenate FASTA
cat $@ > ${OUTP}.fa

# Filter sequences
python ${BASEDIR}/../src/filter.py -d -m ${OUTP}.fa > ${OUTP}.filter.fa
rm ${OUTP}.fa

# Alignment
mafft --auto --thread ${THREADS} --addfragments ${OUTP}.filter.fa ${REF} > ${OUTP}.align
rm ${OUTP}.filter.fa

# Filter alignment
python ${BASEDIR}/../src/filter.py -d -m ${OUTP}.align > ${OUTP}.msa
rm ${OUTP}.align

# Phylogenetic inference
iqtree -nt ${THREADS} -s ${OUTP}.msa
# or: -m GTR+I+G
# or: -B 1000 -alrt 1000

# Plotting or use FigTree
Rscript ${BASEDIR}/../scripts/phylogeny.R ${OUTP}.msa.treefile 
