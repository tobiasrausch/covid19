SHELL := /bin/bash

# Targets
TARGETS = .conda .channels .install .check .pangolin .llama
PBASE=$(shell pwd)

all:   	$(TARGETS)

.conda:
	wget 'https://repo.continuum.io/miniconda/Miniconda3-latest-Linux-x86_64.sh' && bash Miniconda3-latest-Linux-x86_64.sh -b -p ${PBASE}/conda && rm -f Miniconda3-latest-Linux-x86_64.sh && touch .conda

.channels: .conda
	export PATH=${PBASE}/conda/bin:${PATH} && conda config --add channels defaults && conda config --add channels conda-forge && conda config --add channels bioconda && touch .channels

.install: .conda .channels
	export PATH=${PBASE}/conda/bin:${PATH} && conda install -y samtools bcftools bedtools htslib bwa trim-galore fastqc delly kraken2 ivar alfred=0.2.3 seqtk freebayes mafft iqtree cyvcf2 scikit-learn verifybamid && touch .install

.check: .conda .channels .install
	export PATH=${PBASE}/conda/bin:${PATH} && samtools --version && bcftools --version && bedtools --version && bgzip --version && tabix --version && trim_galore --version && delly --version && ivar version && touch .check

.pangolin: .conda .channels .install .check
	export PATH=${PBASE}/conda/bin:${PATH} && git clone --recursive https://github.com/cov-lineages/pangolin.git && cd pangolin && conda env create -f environment.yml && source activate pangolin && python setup.py install && pangolin --update && pangolin -v && pangolin -pv && pangolin -lv && cd ../ && touch .pangolin

.llama: .conda .channels .install .check
	export PATH=${PBASE}/conda/bin:${PATH} && git clone --recursive https://github.com/cov-lineages/llama.git && cd llama && conda env create -f environment.yml && source activate llama && python setup.py install && llama -v && cd ../ && touch .llama

clean:
	rm -rf $(TARGETS) $(TARGETS:=.o) conda/ pangolin/
