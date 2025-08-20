#!/bin/bash

if [[ ! $# -eq 2 ]]; then
  echo "Usage: genes.bootstrap.sh <data_dir> <destination_dir>"
  exit 1
fi
DATA_DIR=$1
INPUT_DIR=$2

aws s3 cp ${DATA_DIR}/bin/magma/NCBI37.3.gene.loc ${INPUT_DIR}/inputs/
