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
THREADS=4

# FastQC
mkdir prefastqc
fastqc -t ${THREADS} -o prefastqc/ ${FQ1}
fastqc -t ${THREADS} -o prefastqc/ ${FQ2}

# Adapter trimming
trim_galore --paired --basename ${OUTP} ${FQ1} ${FQ2} > ${OUTP}.trim_galore.log 2> ${OUTP}.trim_galore.err
