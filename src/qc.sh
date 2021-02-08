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

# Input parameters
FASTA=${1}
OUTP=${2}
REF=${BASEDIR}/../ref/NC_045512.2.fa
THREADS=4

# Alignment to SARS-CoV-2 reference
if [ ! -f ${OUTP}.align.fa.gz ]
then
    # Trim leading and trailing Ns, replace IUPAC ambiguous characters with N
    cat ${FASTA} | sed 's/^N*//' | sed 's/N*$//' | tr 'RYSWKMBDHV' 'NNNNNNNNNN' > ${OUTP}.in.fasta

    # Alignment to the reference (end-gaps free in reference), score matches as 1 for percent identity
    alfred pwalign -q -f h -a ${OUTP}.align.fa.gz -g 0 -e 0 -m 1 -n 0 ${REF} ${OUTP}.in.fasta > ${OUTP}.alistats
    rm ${OUTP}.in.fasta
fi

# Compute summary QC table
python ${BASEDIR}/qc.py -p ${OUTP} > ${OUTP}.qc.summary
