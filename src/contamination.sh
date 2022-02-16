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

source activate covid19

# Input parameters
FQ1=${1}
FQ2=${2}
OUTP=${3}
STDDB=${BASEDIR}/../kraken2/stdDB/
THREADS=4

# Only run if kraken2 standard DB present
if [ -d ${STDDB} ]
then    
    # Kraken2 standard DB filtering
    kraken2 --paired --db ${STDDB} --threads ${THREADS} --output ${OUTP}.kraken2.out.txt --report ${OUTP}.contamination.report.txt ${FQ1} ${FQ2}
    rm ${OUTP}.kraken2.out.txt
fi
