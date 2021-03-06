#!/bin/bash

if [ $# -ne 2 ]
then
    echo ""
    echo "Usage: $0 <input.fasta> <unique_sample_id>"
    echo ""
    exit -1
fi

SCRIPT=$(readlink -f "$0")
BASEDIR=$(dirname "$SCRIPT")
export PATH=${BASEDIR}/../conda/bin:${PATH}

# Input parameters
FASTA=${1}
OUTP=${2}
REF=${BASEDIR}/../ref/NC_045512.2.fa
THREADS=4
RUNVADR=0

# Run pangolin
if [ -f ${FASTA} ]
then
    source activate pangolin
    pangolin -t ${THREADS} --outfile ${OUTP}.lineage.csv ${FASTA}
    conda deactivate

    # Optional: nextclade (requires docker)
    if command -v docker &> /dev/null
    then
	# run nextclade container
	cp ${FASTA} /tmp
	if [ ! -f /tmp/reference.fasta ]
	then
	    wget -O /tmp/reference.fasta https://raw.githubusercontent.com/nextstrain/nextclade/master/data/sars-cov-2/reference.fasta
	fi
	if [ ! -f /tmp/genemap.gff ]
	then
	    wget -O /tmp/genemap.gff https://raw.githubusercontent.com/nextstrain/nextclade/master/data/sars-cov-2/genemap.gff
	fi
	if [ ! -f /tmp/tree.json ]
	then
	    wget -O /tmp/tree.json https://raw.githubusercontent.com/nextstrain/nextclade/master/data/sars-cov-2/tree.json
	fi
	if [ ! -f /tmp/qc.json ]
	then
	    wget -O /tmp/qc.json https://raw.githubusercontent.com/nextstrain/nextclade/master/data/sars-cov-2/qc.json
	fi
	if [ ! "$(docker ps -q -f name=nextclade)" ]; then
	    if [ "$(docker ps -aq -f status=exited -f name=nextclade)" ]; then
		# cleanup
		docker rm nextclade
	    fi
	    # run nextclade container
	    docker run -it --name nextclade --rm -u `id -u` --volume="/tmp/:/seq" nextstrain/nextclade nextclade --input-fasta "/seq/${FASTA}" --input-root-seq "/seq/reference.fasta" --input-tree "/seq/tree.json" --input-gene-map "/seq/genemap.gff" --input-qc-config "/seq/qc.json" --output-basename ${OUTP} --output-dir "/seq/" --output-json /seq/${OUTP}.json --output-tsv /seq/${OUTP}.tsv
	    cp /tmp/${OUTP}.json .
	    cp /tmp/${OUTP}.tsv .
	    rm /tmp/${OUTP}*
	fi

	# VADR: https://github.com/ncbi/vadr
	if [ ${RUNVADR} -eq 1 ]
	then
	    if [ ! "$(docker ps -q -f name=vadr)" ]; then
		if [ "$(docker ps -aq -f status=exited -f name=vadr)" ]; then
		    # cleanup
		    docker rm vadr
		fi
		# run vadr container
		docker run -it --name vadr --rm -u `id -u` --volume="/tmp/:/data" staphb/vadr /opt/vadr/vadr/v-annotate.pl --mxsize 64000 -s -r --nomisc --mkey NC_045512 --lowsim5term 2 --lowsim3term 2 --fstlowthr 0.0 --alt_fail lowscore,fsthicnf,fstlocnf "/data/${FASTA}" ${OUTP}.out
		cp /tmp/${OUTP}.out/${OUTP}.out.vadr.pass.tbl .
		cp /tmp/${OUTP}.out/${OUTP}.out.vadr.fail.tbl .
		rm -rf /tmp/${OUTP}.out/
	    fi
	else
	    touch ${OUTP}.out.vadr.pass.tbl
	    touch ${OUTP}.out.vadr.fail.tbl
	fi
	rm -f /tmp/${FASTA}
    fi
fi
