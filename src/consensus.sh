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

# Align both consensus sequences
alfred pwalign -k -f v -a ${OUTP}.vert.gz -g -4 -e -2 -m 5 -n -4 ${OUTP}.cons.fa ${OUTP}.fb.fa
zcat ${OUTP}.vert.gz | sed 's/^\(.\)/\1 /' | awk '{if ($1!="-") { P1+=1; }; if ($2!="-") {P2+=1;}; print P1"\t"P2"\t"$1"\t"$2"\t"($1!=$2);}' | gzip -c > ${OUTP}.cons.diff.gz
rm ${OUTP}.vert.gz ${OUTP}.fb.fa
