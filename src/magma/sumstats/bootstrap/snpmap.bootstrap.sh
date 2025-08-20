#!/bin/bash

if [[ ! $# -eq 3 ]]; then
  echo "Usage: snpmap.bootstrap.sh <data_dir> <destination_dir> <genome_build>"
  exit 1
fi
DATA_DIR=$1
INPUT_DIR=$2
GENOME_BUILD=$3

aws s3 cp ${DATA_DIR}/bin/magma/snpmap/sumstats.${GENOME_BUILD}.snpmap ${INPUT_DIR}/snpmap/
