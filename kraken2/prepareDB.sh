#!/bin/bash

SCRIPT=$(readlink -f "$0")
BASEDIR=$(dirname "$SCRIPT")
export PATH=${BASEDIR}/../conda/bin:${PATH}

THREADS=12

# Build standard DB, only required for ./src/contamination.sh
#kraken2-build --standard --threads ${THREADS} --db stdDB
#kraken2-build --clean --db stdDB

# Build human DB
mkdir humanDB
kraken2-build --download-taxonomy --db humanDB
kraken2-build --download-library human --threads ${THREADS} --db humanDB
kraken2-build --build --threads ${THREADS} --db humanDB
kraken2-build --clean --db humanDB
