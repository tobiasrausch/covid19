#!/bin/bash

if [ $# -ne 2 ]
then
    echo ""
    echo "Usage: $0 <input.bam> <unique_sample_id>"
    echo ""
    exit -1
fi

SCRIPT=$(readlink -f "$0")
BASEDIR=$(dirname "$SCRIPT")
export PATH=${BASEDIR}/../conda/bin:${PATH}

# Input parameters
BAM=${1}
OUTP=${2}
REF=${BASEDIR}/../ref/NC_045512.2.fa
ANNO=${BASEDIR}/../ref/NC_045512.2.anno.bcf
THREADS=4

# call variants
freebayes --ploidy 1 --no-partial-observations --min-alternate-count 10 --report-genotype-likelihood-max --min-alternate-fraction 0.15 --fasta-reference ${REF} --genotype-qualities -v ${OUTP}.vcf ${BAM}
bgzip ${OUTP}.vcf
tabix ${OUTP}.vcf.gz

# normalize variants
bcftools norm -f ${REF} -O b -o ${OUTP}.norm.bcf ${OUTP}.vcf.gz
rm ${OUTP}.vcf.gz ${OUTP}.vcf.gz.tbi

# filter on quality and depth
bcftools filter -i 'QUAL>20 && INFO/DP>=10' -O b -o ${OUTP}.filtered.bcf ${OUTP}.norm.bcf
bcftools index ${OUTP}.filtered.bcf
rm ${OUTP}.norm.bcf 

# annotate bcf
bcftools annotate -O b -o ${OUTP}.bcf -a ${ANNO} -c INFO/CSQ ${OUTP}.filtered.bcf
bcftools index ${OUTP}.bcf
rm ${OUTP}.filtered.bcf ${OUTP}.filtered.bcf.csi

# compute the depth at each reference position
samtools depth -aa -d 0 ${BAM} > ${OUTP}.depth

# generate tab-delimited variant table (IMPACT >= LOW) and key mutation string (csv)
python ${BASEDIR}/call.py -s ${OUTP} -d ${OUTP}.depth -v ${OUTP}.bcf -o ${OUTP}.mutation.csv > ${OUTP}.variants.tsv
