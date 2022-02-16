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
GRCH38=${BASEDIR}/../ref/Homo_sapiens.GRCh38.dna.primary_assembly.fa
HDB=${BASEDIR}/../kraken2/humanDB/
THREADS=4

# Check human DB
if [ ! -d ${HDB} ]; then echo "${HDB} not found!"; exit; fi

# Map against host reference
bwa mem -R "@RG\tID:${OUTP}\tSM:${OUTP}" -t ${THREADS} ${GRCH38} ${FQ1} ${FQ2} | samtools sort -@ ${THREADS} -o ${OUTP}.srt.bam
samtools index ${OUTP}.srt.bam
samtools view -b ${OUTP}.srt.bam `cut -f 1 ${GRCH38}.fai | grep -v "NC_045512.2" | tr '\n' ' '` > ${OUTP}.grch38.bam
samtools index ${OUTP}.grch38.bam
samtools flagstat ${OUTP}.grch38.bam > ${OUTP}.grch38.bam.flagstat
rm ${OUTP}.srt.bam ${OUTP}.srt.bam.bai

# Filter host reads (mapping to GRCh38)
zcat ${FQ1} | awk 'NR%4==1' | cut -c 2- | sed 's/[ \t].*$//' | sort | uniq > ${OUTP}.all.reads
samtools view ${OUTP}.grch38.bam | cut -f 1 | sort | uniq > ${OUTP}.remove.reads
sort -m ${OUTP}.all.reads ${OUTP}.remove.reads | uniq -u | sed 's/^/@/' > ${OUTP}.retain.reads
wc -l *.reads > ${OUTP}.host.reads
zcat ${FQ1} | grep --no-group-separator -w -Ff ${OUTP}.retain.reads -A 3 | gzip -c > ${OUTP}.grch38.1.fq.gz
zcat ${FQ2} | grep --no-group-separator -w -Ff ${OUTP}.retain.reads -A 3 | gzip -c > ${OUTP}.grch38.2.fq.gz
rm ${OUTP}.all.reads ${OUTP}.remove.reads ${OUTP}.retain.reads
rm ${FQ1} ${FQ2}

# Kraken2 host read filter
kraken2 --paired --db ${HDB} --threads ${THREADS} --output ${OUTP}.kraken2.out.txt --report ${OUTP}.kraken2.report.txt --unclassified-out ${OUTP}.filtered.R#.fq ${OUTP}.grch38.1.fq.gz ${OUTP}.grch38.2.fq.gz
gzip ${OUTP}.filtered.R_1.fq
gzip ${OUTP}.filtered.R_2.fq
rm ${OUTP}.grch38.1.fq.gz ${OUTP}.grch38.2.fq.gz ${OUTP}.kraken2.out.txt

# FastQC
mkdir postfastqc
fastqc -t ${THREADS} -o postfastqc/ ${OUTP}.filtered.R_1.fq.gz
fastqc -t ${THREADS} -o postfastqc/ ${OUTP}.filtered.R_2.fq.gz
