#!/bin/bash

if [[ ! $# -eq 2 ]]; then
  echo "Usage: g1000.bootstrap.sh <data_dir> <destination_dir> <ancestry>"
  exit 1
fi
DATA_DIR=$1
INPUT_DIR=$2

mkdir -p ${INPUT_DIR}/gene_loc
aws s3 cp ${DATA_DIR}/bin/gene_loc/gene.loc ${INPUT_DIR}/gene_loc/
