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

# Trim leading and trailing Ns, replace IUPAC ambiguous characters with N
cat ${FASTA} | sed 's/^N*//' | sed 's/N*$//' | tr 'RYSWKMBDHV' 'NNNNNNNNNN' > ${OUTP}.in.fasta

# Alignment to the reference (end-gaps free in reference)
alfred pwalign -q -f h -a ${OUTP}.align.fa.gz ${REF} ${OUTP}.in.fasta
rm ${OUTP}.in.fasta

# Compute summary QC table
python ${BASEDIR}/qc.py -p ${OUTP} > ${OUTP}.qc.summary
