SHELL := /bin/bash

# Targets
TARGETS = .conda .channels .install .check .pangolin
PBASE=$(shell pwd)

all:   	$(TARGETS)

.conda:
	wget 'https://repo.continuum.io/miniconda/Miniconda3-latest-Linux-x86_64.sh' && bash Miniconda3-latest-Linux-x86_64.sh -b -p ${PBASE}/conda && rm -f Miniconda3-latest-Linux-x86_64.sh && touch .conda

.channels: .conda
	export PATH=${PBASE}/conda/bin:${PATH} && conda config --add channels defaults && conda config --add channels conda-forge && conda config --add channels bioconda && touch .channels

.install: .conda .channels
	export PATH=${PBASE}/conda/bin:${PATH} && conda install -y -c bioconda samtools htslib bcftools python -mn covpipe_env && source activate covpipe_env && pip install covpipe && touch .install

.check: .conda .channels .install
	export PATH=${PBASE}/conda/bin:${PATH} && source activate covpipe_env && which pip && ncov_minipipe -h && touch .check

.pangolin: .conda .channels .install .check
	export PATH=${PBASE}/conda/bin:${PATH} && source activate covpipe_env && pip install matplotlib scikit-learn && update_pangolin && touch .pangolin

clean:
	rm -rf $(TARGETS) $(TARGETS:=.o) conda/
