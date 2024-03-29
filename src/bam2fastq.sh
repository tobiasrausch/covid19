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

source activate	covid19

# Input parameters
BAM=${1}
OUTP=${2}
THREADS=4

# sort by name
samtools sort -@ ${THREADS} -n ${BAM} -o ${OUTP}.tmp.bam

# create paired-end files
samtools fastq -@ ${THREADS} -1 ${OUTP}.R1.fastq.gz -2 ${OUTP}.R2.fastq.gz -0 /dev/null -s ${OUTP}.fastq.gz -n ${OUTP}.tmp.bam
rm ${OUTP}.tmp.bam
