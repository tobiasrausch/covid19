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

source activate	covid19

# Input parameters
FQ1=${1}
FQ2=${2}
OUTP=${3}
REF=${BASEDIR}/../ref/NC_045512.2.fa
PRIMER=${BASEDIR}/../ref/nCoV-2019.primer.bed
THREADS=4

# Alignment
bwa mem -R "@RG\tID:${OUTP}\tSM:${OUTP}" -t ${THREADS} ${REF} ${FQ1} ${FQ2} | samtools sort -@ ${THREADS} -o ${OUTP}.srt.bam
samtools index ${OUTP}.srt.bam
samtools flagstat ${OUTP}.srt.bam > ${OUTP}.srt.bam.flagstat

# Guess amplicon design
BESTCOV=0
cat ${BASEDIR}/../ref/*.primer.bed | cut -f 1-3 | sort | uniq -u > ${OUTP}.fetch
for PRIN in nCoV-2019 neb_vss1a neb_vss2a
do
    if [ -f ${BASEDIR}/../ref/${PRIN}.primer.bed ]
    then
	cat ${BASEDIR}/../ref/${PRIN}.primer.bed | grep -w -Ff ${OUTP}.fetch > ${OUTP}.${PRIN}.unique
	TOTALLEN=`cat ${OUTP}.${PRIN}.unique | awk '{SUM+=$3-$2;} END {print SUM;}'`
	COVPRIN=`samtools depth -a -d 0 -b ${OUTP}.${PRIN}.unique ${OUTP}.srt.bam | awk '{SUM+=$3;} END {print int(SUM/NR);}'`
	echo -e "${PRIN}\t${COVPRIN}\t${TOTALLEN}" >> ${OUTP}.srt.bam.flagstat
	if [ ${COVPRIN} -gt ${BESTCOV} ]
	then
	    PRIMER=${BASEDIR}/../ref/${PRIN}.primer.bed
	    BESTCOV=${COVPRIN}
	fi
	rm ${OUTP}.${PRIN}.unique
    fi
done
rm -f ${OUTP}.fetch
echo "Used" ${PRIMER} >> ${OUTP}.srt.bam.flagstat

# Alignment QC
alfred qc -r ${REF} -j ${OUTP}.alfred.json.gz -o ${OUTP}.alfred.tsv.gz ${OUTP}.srt.bam

# Mask priming regions
ivar trim -e -i ${OUTP}.srt.bam -b ${PRIMER} -p ${OUTP}.trim > ${OUTP}.iVar.trim
samtools sort -@ ${THREADS} -o ${OUTP}.trim.srt.bam ${OUTP}.trim.bam
samtools index ${OUTP}.trim.srt.bam
rm ${OUTP}.trim.bam

# Replace alignments
mv ${OUTP}.trim.srt.bam ${OUTP}.srt.bam
mv ${OUTP}.trim.srt.bam.bai ${OUTP}.srt.bam.bai
