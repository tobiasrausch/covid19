#!/bin/bash

if [ $# -ne 3 ]
then
    echo ""
    echo "Usage: $0 <read.1.fq.gz> <read.2.fq.gz> <unique_sample_id>"
    echo ""
    exit -1
fi

SCRIPT=$(readlink -f "$0")
BASEDIR=$(dirname "$SCRIPT")
export PATH=${BASEDIR}/../conda/bin:${PATH}

# Input parameters
FQ1=${1}
FQ2=${2}
OUTP=${3}
REF=${BASEDIR}/../ref/NC_045512.2.fa
PRIMER=${BASEDIR}/../ref/nCoV-2019.primer.bed
THREADS=4

# Alignment
bwa mem -R "@RG\tID:${OUTP}\tSM:${OUTP}" -t ${THREADS} ${REF} ${FQ1} ${FQ2} | samtools sort -@ ${THREADS} -o ${OUTP}.srt.bam
samtools index ${OUTP}.srt.bam
samtools flagstat ${OUTP}.srt.bam > ${OUTP}.srt.bam.flagstat

# Alignment QC
alfred qc -r ${REF} -j ${OUTP}.alfred.json.gz -o ${OUTP}.alfred.tsv.gz ${OUTP}.srt.bam

# Mask priming regions
ivar trim -e -i ${OUTP}.srt.bam -b ${PRIMER} -p ${OUTP}.trim > ${OUTP}.iVar.trim
samtools sort -@ ${THREADS} -o ${OUTP}.trim.srt.bam ${OUTP}.trim.bam
samtools index ${OUTP}.trim.srt.bam
rm ${OUTP}.trim.bam

# Replace alignments
mv ${OUTP}.trim.srt.bam ${OUTP}.srt.bam
mv ${OUTP}.trim.srt.bam.bai ${OUTP}.srt.bam.bai
