#!/bin/bash

if [[ ! $# -eq 2 ]]; then
  echo "Usage: sldsc.bootstrap.sh <data_dir> <destination_dir>"
  exit 1
fi
DATA_DIR=$1
INPUT_DIR=$2

aws s3 cp ${DATA_DIR}/bin/phenotype_map/phenotype_map.tsv ${INPUT_DIR}/phenotype_map/
