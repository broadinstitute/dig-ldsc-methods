#!/bin/bash

if [[ ! $# -eq 3 ]]; then
  echo "Usage: phenotype_genes.bootstrap.sh <data_dir> <destination_dir> <ancestry>"
  exit 1
fi
DATA_DIR=$1
INPUT_DIR=$2
ANCESTRY=$3

mkdir -p ${INPUT_DIR}/inputs/
aws s3 cp ${DATA_DIR}/bin/magma/magma.${ANCESTRY}.genes.zip ${INPUT_DIR}/inputs/
