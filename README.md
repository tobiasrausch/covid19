# SARS-CoV-2 data analysis

SARS-CoV-2 analysis pipeline for short-read, paired-end sequencing.

## Installation

A `Makefile` is part of the code that installs all dependencies using bioconda.

`make all`

## Preparing the reference databases and indexes

There is a script to download and index the SARS-CoV-2 and GRCh38 reference sequence.

`cd ref/ && ./prepareREF.sh`

There is another script to prepare the kraken2 human database to filter host reads.

`cd kraken2/ && ./prepareDB.sh`

## Running the data analysis pipeline

There is a run script that performs adapter trimming, host read removal, alignment, variant calling and annotation, consensus calling and some quality control (work-in-progress). The last parameter, called `unique_sample_id`, is used to create a unique output prefix and a directory in the current working directory.

`./src/run.sh <read.1.fq.gz> <read.2.fq.gz> <unique_sample_id>`

## Example

The repository contains an example script using a [COG-UK](https://www.cogconsortium.uk/) data set.

`cd example/ && ./expl.sh`

## Credits

Many thanks to the open-science of [COG-UK](https://www.cogconsortium.uk/), their data sets in [ENA](https://www.ebi.ac.uk/ena/browser/home) were very useful to develop the code. The workflow uses many tools distributed via [bioconda](https://bioconda.github.io/), please see the [Makefile](https://github.com/tobiasrausch/covid19/blob/main/Makefile) for all the dependencies and of course, thanks to all the developers.
