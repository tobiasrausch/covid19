SHELL := /bin/bash

# Targets
TARGETS = .check .conda .mamba .tools .pangolin
PBASE=$(shell pwd)

all: ${TARGETS}

.conda:
	wget 'https://repo.anaconda.com/miniconda/Miniconda3-latest-Linux-x86_64.sh' && bash Miniconda3-latest-Linux-x86_64.sh -b -p ${PBASE}/conda && rm -f Miniconda3-latest-Linux-x86_64.sh && touch .conda

.mamba: .conda
	export PATH=${PBASE}/conda/bin:${PATH} && conda install -y -n base -c conda-forge mamba && touch .mamba

.tools: .conda .mamba
	export PATH=${PBASE}/conda/bin:${PATH} && source activate base && mamba create -y -c conda-forge -c bioconda -n covid19 samtools bcftools bedtools htslib bwa trim-galore fastqc delly kraken2 ivar alfred freebayes mafft iqtree cyvcf2 scikit-learn verifybamid wally && touch .tools

.pangolin: .conda .mamba
	export PATH=${PBASE}/conda/bin:${PATH} && wget https://github.com/cov-lineages/pangolin/archive/refs/tags/v4.0.6.tar.gz && tar -xzf v4.0.6.tar.gz && rm v4.0.6.tar.gz && mv pangolin-*/ pangolin && cd pangolin && conda env create -f environment.yml && source activate pangolin && pip install matplotlib scipy && pip install . && pangolin --update && pangolin -v && pangolin -pv && cd ../ && touch .pangolin

.check: .conda .mamba .tools
	export PATH=${PBASE}/conda/bin:${PATH} && source activate covid19 && samtools --version && bcftools --version && bedtools --version && bgzip --version && tabix --version && trim_galore --version && delly --version && ivar version && touch .check

clean:
	rm -rf $(TARGETS) $(TARGETS:=.o) conda/ pangolin/
