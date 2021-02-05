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
export PATH=/opt/dev/covid19/benchmarking/artic-nf/conda/bin:${PATH}

# Parameters
FQ1=${1}   #FASTQ read1
FQ2=${2}   #FASTQ read2
OUTP=${3}  # output prefix
ANNO=${BASEDIR}/../ncov2019-artic-nf/typing/MN908947.3.gff
YML=${BASEDIR}/../ncov2019-artic-nf/typing/SARS-CoV-2.types.yaml 

# Check files
if [ ! -f ${FQ1} ]; then echo "${FQ1} not found!"; exit; fi
if [ ! -f ${FQ2} ]; then echo "${FQ2} not found!"; exit; fi
if [ ! -f ${REF} ]; then echo "${REF} not found!"; exit; fi
if [ ! -f ${ANNO} ]; then echo "${ANNO} not found!"; exit; fi
if [ ! -f ${YML} ]; then echo "${YML} not found!"; exit; fi

# Stage-in the data for scratch, cloud, ...
mkdir -p ${OUTP}
cp ${FQ1} ${OUTP}/${OUTP}_S1_L001_R1_001.fastq.gz
cp ${FQ2} ${OUTP}/${OUTP}_S1_L001_R2_001.fastq.gz
cd ${OUTP}

# artic nf
#nextflow run connor-lab/ncov2019-artic-nf --help
nextflow run connor-lab/ncov2019-artic-nf --illumina --prefix ${OUTP} --directory `pwd` --gff ${ANNO} --allowNoprimer true --illuminaKeepLen 35 --illuminaQualThreshold 20 --ivarFreqThreshold 0.9 --ivarMinDepth 20 --ivarMinFreqThreshold 0.15 --ivarMinVariantQuality 20 --yaml ${YML}
