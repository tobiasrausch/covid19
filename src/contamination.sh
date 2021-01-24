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
STDDB=${BASEDIR}/../kraken2/stdDB/
THREADS=4

# Check standard DB
if [ ! -d ${STDDB} ]; then echo "${STDDB} not found!"; exit; fi

# Kraken2 standard DB filtering
kraken2 --paired --db ${STDDB} --threads ${THREADS} --output ${OUTP}.kraken2.out.txt --report ${OUTP}.kraken2.report.txt --unclassified-out ${OUTP}.filtered.R#.fq ${FQ1} ${FQ2}
gzip ${OUTP}.filtered.R_1.fq
gzip ${OUTP}.filtered.R_2.fq
