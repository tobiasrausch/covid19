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

source activate covid19

# Input parameters
REF=${BASEDIR}/../ref/NC_045512.2.fa
OUTP=${1}
shift 1

# Generate SNP sites with estimated AC/AN
bcftools merge -0 --no-version $@ | bcftools view -m2 -M2 -v snps --min-ac 1 - | cut -f 1-8 | bgzip > ${OUTP}.vcf.gz
tabix ${OUTP}.vcf.gz

# Re-genotype
NUM=1
echo "Genotyping..."
for F in $@
do
    BAM=`echo ${F} | sed 's/.bcf$/.srt.bam/'`
    if [ -f ${BAM} ]
    then
	echo ${BAM}
	freebayes -l -@ ${OUTP}.vcf.gz --ploidy 1 --report-genotype-likelihood-max --fasta-reference ${REF} --genotype-qualities -v ${OUTP}.geno.${NUM}.vcf ${BAM}
	bgzip ${OUTP}.geno.${NUM}.vcf
	tabix ${OUTP}.geno.${NUM}.vcf.gz
	NUM=`expr ${NUM} + 1`
    fi
done
bcftools merge --no-version -O b -o ${OUTP}.bcf ${OUTP}.geno.*.vcf.gz
bcftools index ${OUTP}.bcf
rm ${OUTP}.geno.*.vcf.gz*
rm ${OUTP}.vcf.gz ${OUTP}.vcf.gz.tbi

# Quick contamination assessment
python ${BASEDIR}/crosscontam.py -b ${OUTP}.bcf > ${OUTP}.tsv
rm ${OUTP}.bcf ${OUTP}.bcf.csi
