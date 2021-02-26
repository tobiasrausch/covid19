#!/bin/bash

if [ $# -lt 2 ]
then
    echo ""
    echo "Usage: $0 <outprefix> <sample1.bcf> <sample2.bcf> ... <sampleN.bcf>"
    echo ""
    exit -1
fi

SCRIPT=$(readlink -f "$0")
BASEDIR=$(dirname "$SCRIPT")
export PATH=${BASEDIR}/../conda/bin:${PATH}

# Input parameters
REF=${BASEDIR}/../ref/NC_045512.2.fa
OUTP=${1}
shift 1

# Generate SNP sites with estimated AC/AN
bcftools merge -0 --no-version $@ | bcftools view -m2 -M2 -v snps --min-ac 1 - | cut -f 1-8 | sed 's/^NC_045512.2/chr22/' | bgzip > ${OUTP}.vcf.gz
tabix ${OUTP}.vcf.gz

# Estimate cross-contamination
NUM=1
rm -f ${OUTP}.freemix.tsv
echo "Cross-contamination..."
for F in $@
do
    BAM=`echo ${F} | sed 's/.bcf$/.srt.bam/'`
    ID=`echo ${BAM} | sed 's/^.*\///' | sed 's/.srt.bam$//'`
    if [ -f ${BAM} ]
    then
	samtools view -h ${BAM} | sed 's/NC_045512.2/chr22/' | samtools view -b - > ${OUTP}.tmp.${NUM}.bam
	samtools index ${OUTP}.tmp.${NUM}.bam
	verifyBamID --ignoreRG --chip-none --site --vcf ${OUTP}.vcf.gz --bam ${OUTP}.tmp.${NUM}.bam --noPhoneHome --maxDepth 1000 --precise --out ${OUTP}.contam.${NUM}
	rm ${OUTP}.tmp.${NUM}.bam ${OUTP}.tmp.${NUM}.bam.bai
	FREEMIX=`cut -f 7 ${OUTP}.contam.${NUM}.selfSM | tail -n 1`
	echo -e "${ID}\t${FREEMIX}" >> ${OUTP}.freemix.tsv
	rm ${OUTP}.contam.${NUM}.*
	NUM=`expr ${NUM} + 1`
    fi
done
