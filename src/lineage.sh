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

source activate covid19

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
    pangolin -t ${THREADS} --analysis-mode accurate --outfile ${OUTP}.lineage.csv ${FASTA}
    conda deactivate

    # Optional: nextclade (requires docker)
    if command -v docker &> /dev/null
    then
	# run nextclade container
	cp ${FASTA} /tmp
	if [ ! -d /tmp/sars-cov-2 ]
	then
	    docker run -it --name nextclade --rm -u `id -u` --volume="/tmp/:/seq" nextstrain/nextclade nextclade dataset get --name sars-cov-2 --output-dir '/seq/sars-cov-2'
	fi
	if [ "$(docker ps -q -f name=nextclade)" ]
	then
	    docker rm -vf nextclade
	    sleep 1
	fi
	docker run -it --name nextclade --rm -u `id -u` --volume="/tmp/:/seq" nextstrain/nextclade nextclade run --input-dataset '/seq/sars-cov-2' --output-basename ${OUTP} --output-all "/seq/" --output-json /seq/${OUTP}.json --output-tsv /seq/${OUTP}.tsv "/seq/${FASTA}"
	cp /tmp/${OUTP}.json .
	cp /tmp/${OUTP}.tsv .
	rm /tmp/${OUTP}*

	# VADR: https://github.com/ncbi/vadr
	if [ ${RUNVADR} -eq 1 ]
	then
	    if [ "$(docker ps -q -f name=vadr)" ]
	    then
		docker rm -vf vadr
	    fi
	    # run vadr container
	    docker run -it --name vadr --rm -u `id -u` --volume="/tmp/:/data" staphb/vadr /opt/vadr/vadr/v-annotate.pl --mxsize 64000 -s -r --nomisc --mkey NC_045512 --lowsim5term 2 --lowsim3term 2 --fstlowthr 0.0 --alt_fail lowscore,fsthicnf,fstlocnf "/data/${FASTA}" ${OUTP}.out
	    cp /tmp/${OUTP}.out/${OUTP}.out.vadr.pass.tbl .
	    cp /tmp/${OUTP}.out/${OUTP}.out.vadr.fail.tbl .
	    rm -rf /tmp/${OUTP}.out/
	else
	    touch ${OUTP}.out.vadr.pass.tbl
	    touch ${OUTP}.out.vadr.fail.tbl
	fi
	rm -f /tmp/${FASTA}
    fi
fi
