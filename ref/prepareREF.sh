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
