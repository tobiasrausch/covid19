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
	docker run -it --name nextclade --rm -u `id -u` --volume="/tmp/:/seq" nextstrain/nextclade nextclade run --input-fasta "/seq/${FASTA}" --input-dataset '/seq/sars-cov-2' --output-basename ${OUTP} --output-dir "/seq/" --output-json /seq/${OUTP}.json --output-tree /seq/${OUTP}.tree.json --output-tsv /seq/${OUTP}.tsv
	cp /tmp/${OUTP}.json .
	cp /tmp/${OUTP}.tsv .
	cp /tmp/${OUTP}.tree.json .
	rm /tmp/${OUTP}*
    fi
fi
