#!/bin/bash

if [ $# -lt 2 ]
then
    echo ""
    echo "Usage: $0 <outprefix> <sample1.depth> <sample2.depth> ... <sampleN.depth>"
    echo ""
    exit -1
fi

SCRIPT=$(readlink -f "$0")
BASEDIR=$(dirname "$SCRIPT")
export PATH=${BASEDIR}/../conda/bin:${PATH}

# Input parameters
OUTP=${1}
shift 1

# Initialize
for F in $@
do
    cat ${F} | cut -f 2 > ${OUTP}.depth
    break
done

# Collect 0x positions
for F in $@
do
    QC=`echo ${F} | sed 's/.depth$/.qc.summary/'`
    if [ -f ${QC} ]
    then
	if [ `grep -c "outcome pass" ${QC}` -eq 1 ]
	then
	    # Stratify further if needed
	    VAR=`echo ${F} | sed 's/.depth$/.variants.tsv/'`
	    #if [ `grep -c "Lineage B.1.351" ${QC}` -eq 1 ]
	    if [ `grep -c -P "22879\tC\tA" ${VAR}` -eq 1 ]
	    then
		echo ${F}
		# 0x
		#cat ${F} | awk '$3==0' | cut -f 2 >> ${OUTP}.depth
		# <10x
		cat ${F} | awk '$3<10' | cut -f 2 >> ${OUTP}.depth
	    fi
	fi
    fi
done

# Aggregate 0x or <10x positions
cat ${OUTP}.depth  | sort | uniq -c | sed 's/^[ \t]*//' | awk '{print $2"\t"$1;}' | sort -n > ${OUTP}.depth.tmp
mv ${OUTP}.depth.tmp ${OUTP}.depth

# Plot
Rscript ${BASEDIR}/../scripts/depth.R ${OUTP}.depth
