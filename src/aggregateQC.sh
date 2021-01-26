#!/bin/bash

if [ $# -lt 2 ]
then
    echo ""
    echo "Usage: $0 <outprefix> <sample1.qc.summary> <sample2.qc.summary> ... <sampleN.qc.summary>"
    echo ""
    exit -1
fi

SCRIPT=$(readlink -f "$0")
BASEDIR=$(dirname "$SCRIPT")
export PATH=${BASEDIR}/../conda/bin:${PATH}

# Input parameters
OUTP=${1}
shift 1

rm -f ${OUTP}.aggr.qc
for F in $@
do
    cat ${F} | datamash -t ' ' transpose >> ${OUTP}.aggr.qc
done
cat ${OUTP}.aggr.qc | sort -r | uniq > ${OUTP}.aggr.qc.tmp
mv ${OUTP}.aggr.qc.tmp ${OUTP}.aggr.qc
