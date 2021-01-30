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
THREADS=4

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
	LP=`pwd`
	if [ ! "$(docker ps -q -f name=nextclade)" ]; then
	    if [ "$(docker ps -aq -f status=exited -f name=nextclade)" ]; then
		# cleanup
		docker rm nextclade
	    fi
	    # run your container
	    docker run -it --name nextclade --rm -u `id -u` --volume="${LP}/:/seq" neherlab/nextclade nextclade --input-fasta "/seq/${FASTA}" --output-json "/seq/${OUTP}.json"
	fi
    fi
fi
