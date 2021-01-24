SHELL := /bin/bash

# Targets
TARGETS = .conda .channels .install .check
PBASE=$(shell pwd)

all:   	$(TARGETS)

.conda:
	wget 'https://repo.continuum.io/miniconda/Miniconda3-latest-Linux-x86_64.sh' && bash Miniconda3-latest-Linux-x86_64.sh -b -p ${PBASE}/conda && rm -f Miniconda3-latest-Linux-x86_64.sh && touch .conda

.channels: .conda
	export PATH=${PBASE}/conda/bin:${PATH} && conda config --add channels defaults && conda config --add channels conda-forge && conda config --add channels bioconda && touch .channels

.install: .conda .channels
	export PATH=${PBASE}/conda/bin:${PATH} && conda install -y samtools bcftools bedtools htslib bwa trim-galore fastqc delly kraken2 biobambam ivar alfred seqtk freebayes mafft cyvcf2 && touch .install

.check: .conda .channels .install
	export PATH=${PBASE}/conda/bin:${PATH} && samtools --version && bcftools --version && bedtools --version && bgzip --version && tabix --version && trim_galore --version && delly --version && ivar version && touch .check

clean:
	rm -rf $(TARGETS) $(TARGETS:=.o) bin/
