# SARS-CoV-2 data analysis

SARS-CoV-2 analysis pipeline for short-read, paired-end sequencing.

## Installation

A [Makefile](https://github.com/tobiasrausch/covid19/blob/main/Makefile) is part of the code that installs all dependencies using bioconda.

`git clone --recursive https://github.com/tobiasrausch/covid19.git`

`cd covid19`

`make all`

## Preparing the reference databases and indexes

There is a script to download and index the SARS-CoV-2 and GRCh38 reference sequence.

`cd ref/ && ./prepareREF.sh`

There is another script to prepare the kraken2 human database to filter host reads.

`cd kraken2/ && ./prepareDB.sh`

## Running the data analysis pipeline

There is a run script that performs adapter trimming, host read removal, alignment, variant calling and annotation, consensus calling and some quality control. The last parameter, called `unique_sample_id`, is used to create a unique output directory in the current working directory.

`./src/run.sh <read.1.fq.gz> <read.2.fq.gz> <unique_sample_id>`

## Output

The main output files are:

* The adapter-trimmed and host-filtered FASTQ files: `ls <unique_sample_id>.filtered.R_[12].fq.gz`

* The alignment to SARS-CoV-2: `ls <unique_sample_id>.srt.bam`

* The consensus sequence: `ls <unique_sample_id>.cons.fa`

* The annotated variants: `ls <unique_sample_id>.variants.tsv`

* The assigned lineage: `ls <unique_sample_id>.lineage.csv`

* The summary QC report: `ls <unique_sample_id>.qc.summary`

## Aggregating results

The above pipeline generates a report for every sample. It can be naively parallelized on the sample level. You can then aggregate all the QC information and the lineage & clade assignments using

`./src/aggregate.sh outtable */*.qc.summary`

## Estimating cross-contamination

You can estimate cross-contamination based on the allelic frequencies of variant calls using

`./src/crosscontam.sh contam */*.bcf`

This works best on good quality consensus sequences, i.e.:

`./src/crosscontam.sh contam `grep "RKI pass" */*.qc.summary | sed 's/.qc.summary.*$/.bcf/' | tr '\n' ' '`

## Example

The repository contains an example script using a [COG-UK](https://www.cogconsortium.uk/) data set.

`cd example/ && ./expl.sh`

## Credits

Many thanks to the open-science of [COG-UK](https://www.cogconsortium.uk/), their data sets in [ENA](https://www.ebi.ac.uk/ena/browser/home) were very useful to develop the code. The workflow uses many tools distributed via [bioconda](https://bioconda.github.io/), please see the [Makefile](https://github.com/tobiasrausch/covid19/blob/main/Makefile) for all the dependencies and of course, thanks to all the developers.
