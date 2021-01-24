# SARS-CoV-2 data analysis

SARS-CoV-2 analysis pipeline for short-read, paired-end illumina sequencing

## Installation

A `Makefile` is part of the code that installs all dependencies using bioconda.

`make all`

## Preparing the reference databases and indexes

There is a script to download and index the SARS-CoV-2 and GRCh38 reference sequence.

`cd ref/ && ./prepareREF.sh`

There is another script to prepare the kraken2 human database to filter host reads.

`cd kraken2/ && ./prepareDB.sh`
