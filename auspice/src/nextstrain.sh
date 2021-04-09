#!/bin/bash

if [ $# -lt 2 ]
then
    echo ""
    echo "Usage: $0 <global|europe> <outprefix> <sample1.cons.fa> <sample2.cons.fa> ... <sampleN.cons.fa>"
    echo ""
    exit -1
fi

SCRIPT=$(readlink -f "$0")
BASEDIR=$(dirname "$SCRIPT")
export PATH=${BASEDIR}/../conda/bin:${PATH}

# Input parameters
REGION=${1}
shift 1
OUTP=${1}
shift 1
DATA=${BASEDIR}/../../../covid/nextstrain/
THREADS=4

# Prepare data
cp ${DATA}/ncov_${REGION}.fasta ${BASEDIR}/../ncov/data/example_sequences.fasta
cp ${DATA}/ncov_${REGION}.tsv ${BASEDIR}/../ncov/data/example_metadata.tsv
for F in $@
do
    echo "Processing..." ${F}
    ID=`echo ${F} | sed 's/^.*\///' | sed 's/.cons.fa//'`
    FOLDER=`echo ${F} | sed 's/\/[^\/]*$//'`
    if [ -f ${FOLDER}/${ID}.qc.summary ]
    then
	OUTSTR=""
	LEN=`cat ${F} | grep -v "^>" | tr -d '\n' | awk '{print length($1);}'`
	CLADE=`cat ${FOLDER}/${ID}.qc.summary | grep -w "^Clade" | cut -f 2 -d ' '`
	LINEAGE=`cat ${FOLDER}/${ID}.qc.summary | grep -w "^Lineage" | cut -f 2 -d ' '`
	DATE=`date +%Y-%m-%d`
	for KEY in `head -n 1 ${BASEDIR}/../ncov/data/example_metadata.tsv`
	do
	    if [ ${KEY} == "strain" ]
	    then
		OUTSTR=${ID}
	    elif [ ${KEY} == "virus" ]
	    then
		OUTSTR=${OUTSTR}"\tbetacoronavirus"
	    elif [ ${KEY} == "region" ]
	    then
		OUTSTR=${OUTSTR}"\tEurope"
	    elif [ ${KEY} == "country" ]
	    then
		OUTSTR=${OUTSTR}"\tGermany"
	    elif [ ${KEY} == "division" ]
	    then
		OUTSTR=${OUTSTR}"\tRhineNeckar"
	    elif [ ${KEY} == "segment" ]
	    then
		OUTSTR=${OUTSTR}"\tgenome"
	    elif [ ${KEY} == "length" ]
	    then
		OUTSTR=${OUTSTR}"\t"${LEN}
	    elif [ ${KEY} == "host" ]
	    then
		OUTSTR=${OUTSTR}"\tHuman"
	    elif [ ${KEY} == "Nextstrain_clade" ]
	    then
		OUTSTR=${OUTSTR}"\t"${CLADE}
	    elif [ ${KEY} == "pango_lineage" ]
	    then
		OUTSTR=${OUTSTR}"\t"${LINEAGE}
	    elif [ ${KEY} == "originating_lab" ]
	    then
		OUTSTR=${OUTSTR}"\tEMBLGeneCore"
	    elif [ ${KEY} == "submitting_lab" ]
	    then
		OUTSTR=${OUTSTR}"\tEMBLGeneCore"
	    elif [ ${KEY} == "date_submitted" ]
	    then
		OUTSTR=${OUTSTR}"\t"${DATE}
	    else
		OUTSTR=${OUTSTR}"\t?"
	    fi
	done
	#echo -e ${OUTSTR} >> ${BASEDIR}/../ncov/data/example_metadata.tsv
	#cat ${F} >> ${BASEDIR}/../ncov/data/example_sequences.fasta
    fi
done

# Run ncov
cd ${BASEDIR}/../ncov/
rm -rf auspice benchmarks logs results
export AUGUR_RECURSION_LIMIT=10000
snakemake --cores 4 --profile ./my_profiles/getting_started

