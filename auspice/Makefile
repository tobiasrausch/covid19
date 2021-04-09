SHELL := /bin/bash

# Targets
TARGETS = .conda .channels .install .check .ncov
PBASE=$(shell pwd)

all:   	$(TARGETS)

.conda:
	wget 'https://repo.anaconda.com/miniconda/Miniconda3-latest-Linux-x86_64.sh' && bash Miniconda3-latest-Linux-x86_64.sh -b -p ${PBASE}/conda && rm -f Miniconda3-latest-Linux-x86_64.sh && touch .conda

.channels: .conda
	export PATH=${PBASE}/conda/bin:${PATH} && conda config --add channels defaults && conda config --add channels conda-forge && conda config --add channels bioconda && touch .channels

.install: .conda .channels
	export PATH=${PBASE}/conda/bin:${PATH} && conda install -y augur nextalign nextstrain-cli nodejs git snakemake && pip install awscli rethinkdb==2.3.0.post6 && npm install --global auspice && nextstrain check-setup --set-default && touch .install

.check: .conda .channels .install
	export PATH=${PBASE}/conda/bin:${PATH} && augur -h && auspice -h && nextstrain -h && touch .check

.ncov: .conda .channels .install .check
	export PATH=${PBASE}/conda/bin:${PATH} && git clone https://github.com/nextstrain/ncov.git && cd ncov && zcat data/example_sequences.fasta.gz > data/example_sequences.fasta && snakemake --cores 4 --profile ./my_profiles/getting_started && touch .ncov

clean:
	rm -rf $(TARGETS) $(TARGETS:=.o) conda/ ncov/
