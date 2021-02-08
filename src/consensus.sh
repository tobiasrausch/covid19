#!/bin/bash

if [ $# -ne 2 ]
then
    echo ""
    echo "Usage: $0 <input.bam> <unique_sample_id>"
    echo ""
    exit -1
fi

SCRIPT=$(readlink -f "$0")
BASEDIR=$(dirname "$SCRIPT")
export PATH=${BASEDIR}/../conda/bin:${PATH}

# Input parameters
BAM=${1}
OUTP=${2}
REF=${BASEDIR}/../ref/NC_045512.2.fa
THREADS=4

# Consensus (strict)
samtools mpileup -aa -A -d 0 -B -Q 0 ${BAM} | ivar consensus -t 0.9 -m 20 -n N -p ${OUTP}.cons

# Fix header
sed -i "s/^>.*$/>${OUTP}/" ${OUTP}.cons.fa

# Consensus computation using freebayes variants (lenient)
cat ${REF} | bcftools consensus ${OUTP}.bcf | sed -e "s/^>.*$/>${OUTP}/" > ${OUTP}.fb.fa

# Replace ambiguous codes, align both consensus sequences and compute diff (ignoring Ns)
cat ${OUTP}.cons.fa | tr 'RYSWKMBDHV' 'NNNNNNNNNN' > ${OUTP}.in.fasta
alfred pwalign -f v -a ${OUTP}.align.fa.gz -g -4 -e -2 -m 5 -n -4 ${OUTP}.in.fasta ${OUTP}.fb.fa
zcat ${OUTP}.align.fa.gz | sed 's/^\(.\)/\1 /' | awk '{print NR"\t"$1"\t"$2"\t"($1!=$2);}' | grep -P "[ACGT\-]\t[ACGT\-]\t1$" -C 5 > ${OUTP}.cons.diff
rm ${OUTP}.align.fa.gz ${OUTP}.in.fasta
