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
source activate covpipe_env

# Parameters
FQ1=${1}   #FASTQ read1
FQ2=${2}   #FASTQ read2
OUTP=${3}  # output prefix
REF=${BASEDIR}/../../../ref/NC_045512.2.fa
ANNO=${BASEDIR}/../../../ref/GCF_009858895.2_ASM985889v3_genomic.gff
KRAKEN2=${BASEDIR}/../kraken2/GRCh28.p13_GBcovid19-2020-05-22/
PANGOLIN=`ls -d ${BASEDIR}/../conda/envs/covpipe_env/lib/*/site-packages/covpipe/data/pangolin`
PRIMER=${BASEDIR}/../primer.bedpe

# Check files
if [ ! -f ${FQ1} ]; then echo "${FQ1} not found!"; exit; fi
if [ ! -f ${FQ2} ]; then echo "${FQ2} not found!"; exit; fi
if [ ! -f ${REF} ]; then echo "${REF} not found!"; exit; fi
if [ ! -f ${ANNO} ]; then echo "${ANNO} not found!"; exit; fi
if [ ! -d ${KRAKEN2} ]; then echo "${KRAKEN2} not found!"; exit; fi
if [ ! -d ${PANGOLIN} ]; then echo "${PANGOLIN} not found!"; exit; fi

# Stage-in the data for scratch, cloud, ...
mkdir -p ${OUTP}
cp ${FQ1} ${OUTP}/${OUTP}.1.fq.gz
cp ${FQ2} ${OUTP}/${OUTP}.2.fq.gz
cd ${OUTP}

# Run ncov
ncov_minipipe --adapter Nextera --annotation_gff ${ANNO} --pangolin ${PANGOLIN} --kraken ${KRAKEN2} --reference ${REF} --primer ${PRIMER} --input `pwd` -o `pwd`
