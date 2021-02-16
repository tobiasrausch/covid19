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
    # Compute alignment to reference
    alfred pwalign -k -f v -a ${OUTP}.align.gz -g -4 -e -2 -m 5 -n -4 ${REF} ${FASTA}
    zcat ${OUTP}.align.gz | sed 's/^\(.\)/\1 /' | awk '{if ($1!="-") { P1+=1; }; if ($2!="-") {P2+=1;}; print P1"\t"P2"\t"$1"\t"$2"\t"($1!=$2);}' | gzip -c > ${OUTP}.align.fa.gz
    rm ${OUTP}.align.gz
fi

# Compute summary QC table
python ${BASEDIR}/qc.py -p ${OUTP} > ${OUTP}.qc.summary
