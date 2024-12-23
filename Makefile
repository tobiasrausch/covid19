SHELL := /bin/bash

# Targets
TARGETS = .mamba .tools .pangolin .check .update
PBASE=$(shell pwd)

all: ${TARGETS}

.mamba:
	curl -L -O "https://github.com/conda-forge/miniforge/releases/latest/download/Miniforge3-$(shell uname)-$(shell uname -m).sh" && bash Miniforge3-$(shell uname)-$(shell uname -m).sh -b -p conda && rm "Miniforge3-$(shell uname)-$(shell uname -m).sh" && touch .mamba

.tools: .mamba
	export PATH=${PBASE}/conda/bin:${PATH} && mamba create -y --override-channels -c conda-forge -c bioconda -n covid19 samtools bcftools bedtools htslib bwa trim-galore fastqc delly kraken2 ivar alfred freebayes mafft iqtree cyvcf2 scikit-learn verifybamid wally && touch .tools

.pangolin: .mamba
	export PATH=${PBASE}/conda/bin:${PATH} && wget https://github.com/cov-lineages/pangolin/archive/refs/tags/v4.3.1.tar.gz && tar -xzf v4.3.1.tar.gz && rm v4.3.1.tar.gz && mv pangolin-*/ pangolin && cd pangolin && cat environment.yml | grep -v "defaults" > environment.yml.tmp && mv environment.yml.tmp environment.yml && mamba env create -f environment.yml && source activate pangolin && cd ../ && touch .pangolin

.update: .mamba .pangolin
	export PATH=${PBASE}/conda/bin:${PATH} && source activate pangolin && cd pangolin/ && pip install . && pangolin --update && pangolin -v && pangolin -pv && cd ../ && touch .update

.check: .mamba .tools
	export PATH=${PBASE}/conda/bin:${PATH} && source activate covid19 && samtools --version && bcftools --version && bedtools --version && bgzip --version && tabix --version && trim_galore --version && delly --version && ivar version && touch .check

clean:
	rm -rf $(TARGETS) $(TARGETS:=.o) conda/ pangolin/
