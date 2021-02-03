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
GFF=${BASEDIR}/../ref/GCF_009858895.2_ASM985889v3_genomic.gff
THREADS=4

# Call variants
samtools mpileup -aa -A -d 10000 -B -Q 0 ${BAM} | ivar variants -r ${REF} -g ${GFF} -m 10 -q 20 -t 0.15 -p ${OUTP}.iVar

# Consensus
samtools mpileup -aa -A -d 10000 -B -Q 0 ${BAM} | ivar consensus -t 0.7 -m 10 -n N -p ${OUTP}.cons

# Fix header
sed -i "s/^>.*$/>${OUTP}/" ${OUTP}.cons.fa

# Nucleotide composition
tail -n +2 ${OUTP}.cons.fa | sed 's/\(.\)/\1\n/g' | grep "." | sort | uniq -c > ${OUTP}.cons.comp

# consensus computation using freebayes variants
cat ${REF} | bcftools consensus ${OUTP}.bcf | sed -e "s/^>.*$/>${OUTP}/" > ${OUTP}.fbvar.fa
bwa index ${OUTP}.fbvar.fa
bwa mem -R "@RG\tID:${OUTP}\tSM:${OUTP}" -t ${THREADS} ${OUTP}.fbvar.fa ${OUTP}.filtered.R_1.fq.gz ${OUTP}.filtered.R_2.fq.gz | samtools sort -@ ${THREADS} -o ${OUTP}.fb.bam
samtools index ${OUTP}.fb.bam

# mask low coverage
head -n 1 ${OUTP}.cons.fa > ${OUTP}.fb.fa
paste <(tail -n +2 ${OUTP}.fbvar.fa | sed 's/\(.\)/\1\n/g' | grep ".") <(samtools depth -aa -d 0 ${OUTP}.fb.bam) | awk '{if ($4<10) {print "N"} else {print $1;} }'  | tr '\n' '#' | sed -e 's/#//g' >> ${OUTP}.fb.fa 
rm ${OUTP}.fbvar.fa* ${OUTP}.fb.bam ${OUTP}.fb.bam.bai

# Compute diff of FreeBayes to iVar consensus
diff <(tail -n +2 ${OUTP}.cons.fa | sed 's/\(.\)/\1\n/g' | grep ".") <(tail -n +2 ${OUTP}.fb.fa | sed 's/\(.\)/\1\n/g' | grep ".") > ${OUTP}.cons.diff
