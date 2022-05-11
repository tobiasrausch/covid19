#!/bin/bash

if [ $# -lt 2 ]
then
    echo ""
    echo "Usage: $0 <outprefix> <sample1.depth> <sample2.depthf> ... <sampleN.depth>"
    echo ""
    exit -1
fi

SCRIPT=$(readlink -f "$0")
BASEDIR=$(dirname "$SCRIPT")
export PATH=${BASEDIR}/../conda/bin:${PATH}

source activate covid19

# Input parameters
REF=${BASEDIR}/../ref/NC_045512.2.fa
OUTP=${1}
shift 1

# Concatenate
paste $@ > ${OUTP}.depth
COLS=`cat ${OUTP}.depth | awk '{print NF;}' | head -n 1`
COLS=`seq 3 3 ${COLS} | tr '\n' ',' | sed 's/,$//'`
cut -f 1,2,${COLS} ${OUTP}.depth > ${OUTP}.depth.tmp
mv ${OUTP}.depth.tmp ${OUTP}.depth
Rscript ../scripts/cumcov.R ${OUTP}.depth

source deactivate
