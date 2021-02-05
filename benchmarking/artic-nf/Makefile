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
	export PATH=${PBASE}/conda/bin:${PATH} && conda install -y python=3 biopython=1.74 libxcb matplotlib=3.3.3 pandas=0.23.0=py36_1 bwa=0.7.17=pl5.22.0_2 samtools=1.10 bcftools=1.10 trim-galore=0.6.5 ivar=1.3 pyvcf=0.6.8 pyyaml=5.3.1 nextflow git && touch .install

.check: .conda .channels .install
	export PATH=${PBASE}/conda/bin:${PATH} && nextflow help && git clone --recursive https://github.com/connor-lab/ncov2019-artic-nf.git && git clone --recursive https://github.com/artic-network/artic-ncov2019.git && touch .check

clean:
	rm -rf $(TARGETS) $(TARGETS:=.o) conda/ ncov2019-artic-nf/ artic-ncov2019/
