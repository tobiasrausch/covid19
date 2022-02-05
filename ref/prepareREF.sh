#!/bin/bash

SCRIPT=$(readlink -f "$0")
BASEDIR=$(dirname "$SCRIPT")
export PATH=${BASEDIR}/../conda/bin:${PATH}

# Index SARS-CoV-2
samtools faidx NC_045512.2.fa
bwa index NC_045512.2.fa

# Download GRCh38+SARS-CoV-2 and index
wget ftp://ftp.ensembl.org/pub/release-102/fasta/homo_sapiens/dna/Homo_sapiens.GRCh38.dna.primary_assembly.fa.gz
gunzip Homo_sapiens.GRCh38.dna.primary_assembly.fa.gz
cat NC_045512.2.fa >> Homo_sapiens.GRCh38.dna.primary_assembly.fa
samtools faidx Homo_sapiens.GRCh38.dna.primary_assembly.fa
bwa index Homo_sapiens.GRCh38.dna.primary_assembly.fa

# Precompute VEP for all single-nucleotide variants (already done)
# Annotation resources: https://covid-19.ensembl.org/info/data/ftp/index.html
# bcftools view NC_045512.2.anno.bcf | grep -v "^#" | cut -f 1-5 > ~/vep/input/input.vcf 
# docker run -t -i -v /home/rausch/vep:/opt/vep/.vep ensemblorg/ensembl-vep ./vep --species sars_cov_2 --cache_version 101 --check_existing --transcript_version --no_stats --cache --offline --compress_output bgzip --format vcf --vcf --force_overwrite --dir_cache /opt/vep/.vep/ --dir_plugins /opt/vep/.vep/Plugins/ --input_file /opt/vep/.vep/input/input.vcf --output_file /opt/vep/.vep/output/out.vcf.gz
