#!/bin/bash

if [ $# -lt 2 ]
then
    echo ""
    echo "Usage: $0 <outprefix> <sample1.cons.fa> <sample2.cons.fa> ... <sampleN.cons.fa>"
    echo ""
    exit -1
fi

SCRIPT=$(readlink -f "$0")
BASEDIR=$(dirname "$SCRIPT")
export PATH=${BASEDIR}/../conda/bin:${PATH}

# Input parameters
OUTP=${1}
shift 1
REF=${BASEDIR}/../ref/NC_045512.2.fa
DATA=${BASEDIR}/../../covid/cog_uk/data/
THREADS=4

# Prepare data
cat $@ > ${OUTP}.fa
echo "name,Lineage,Clade" > ${OUTP}.csv
for F in $@
do
    QC=`echo ${F} | sed 's/.cons.fa/.qc.summary/'`
    cat ${QC} | sed 's/ /\t/' | datamash transpose | cut -f 1,4,5 | sed 's/\t/,/g' | tail -n 1 >> ${OUTP}.csv
done

# Run llama
source activate llama
llama -r -t ${THREADS} --colour-fields Clade --label-fields Lineage -i ${OUTP}.csv -f ${OUTP}.fa -d ${DATA} -o ${OUTP}
conda deactivate
