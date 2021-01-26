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

# Parameters
FQ1=${1}   #FASTQ read1
FQ2=${2}   #FASTQ read2
OUTP=${3}  # output prefix
REF=${BASEDIR}/../ref/NC_045512.2.fa
ANNO=${BASEDIR}/../ref/GCF_009858895.2_ASM985889v3_genomic.gff

# Check files
if [ ! -f ${FQ1} ]; then echo "${FQ1} not found!"; exit; fi
if [ ! -f ${FQ2} ]; then echo "${FQ2} not found!"; exit; fi
if [ ! -f ${REF} ]; then echo "${REF} not found!"; exit; fi
if [ ! -f ${ANNO} ]; then echo "${ANNO} not found!"; exit; fi

# Stage-in the data for scratch, cloud, ...
mkdir -p ${OUTP}
cp ${FQ1} ${OUTP}/${OUTP}.1.fq.gz
cp ${FQ2} ${OUTP}/${OUTP}.2.fq.gz
cd ${OUTP}

# Trim adapters
${BASEDIR}/trim.sh ${OUTP}.1.fq.gz ${OUTP}.2.fq.gz ${OUTP}

# Unbiased contamination assessment (optional: requires kraken2 standard DB)
${BASEDIR}/contamination.sh ${OUTP}_val_1.fq.gz ${OUTP}_val_2.fq.gz ${OUTP}

# Clean host reads
${BASEDIR}/clean.sh ${OUTP}_val_1.fq.gz ${OUTP}_val_2.fq.gz ${OUTP}

# Alignment
${BASEDIR}/align.sh ${OUTP}.filtered.R_1.fq.gz ${OUTP}.filtered.R_2.fq.gz ${OUTP}

# Variant calling
${BASEDIR}/call.sh ${OUTP}.srt.bam ${OUTP}

# Consensus computation
${BASEDIR}/consensus.sh ${OUTP}.srt.bam ${OUTP}

# Quality control
${BASEDIR}/qc.sh ${OUTP}.cons.fa ${OUTP}

# Clean-up
rm ${OUTP}.1.fq.gz ${OUTP}.2.fq.gz
