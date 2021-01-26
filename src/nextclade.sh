#!/bin/bash

if [ $# -lt 2 ]
then
    echo ""
    echo "Usage: $0 <outprefix> <sample1.cons.fa> <sample2.cons.fa> ... <sampleN.cons.fa>"
    echo ""
    exit -1
fi

SCRIPT=$(readlink -f "$0")
BASEDIR=$(dirname "$SCRIPT")
export PATH=${BASEDIR}/../conda/bin:${PATH}

# Input parameters
OUTP=${1}
shift 1
THREADS=4

# check presence of docker
if ! command -v docker &> /dev/null
then
    echo "docker could not be found!"
    exit;
fi

# run nextclade container
LP=`pwd`
cat $@ > ${OUTP}.fasta
if [ ! "$(docker ps -q -f name=nextclade)" ]; then
    if [ "$(docker ps -aq -f status=exited -f name=nextclade)" ]; then
        # cleanup
        docker rm nextclade
    fi
    # run your container
    docker run -it --name nextclade --rm -u `id -u` --volume="${LP}/:/seq" neherlab/nextclade nextclade --input-fasta "/seq/${OUTP}.fasta" --output-json "/seq/${OUTP}.json"
fi
