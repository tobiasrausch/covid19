#!/bin/bash

# Download FASTQs
if [ ! -f ERR4867056_1.fastq.gz ]
then
    wget ftp://ftp.sra.ebi.ac.uk/vol1/fastq/ERR486/006/ERR4867056/ERR4867056_1.fastq.gz
fi
if [ ! -f ERR4867056_2.fastq.gz ]
then
    wget ftp://ftp.sra.ebi.ac.uk/vol1/fastq/ERR486/006/ERR4867056/ERR4867056_2.fastq.gz
fi

# Run the analysis pipeline
if [ ! -d sampleXYZ ]
then
    ../src/run.sh ERR4867056_1.fastq.gz ERR4867056_2.fastq.gz sampleXYZ
fi
