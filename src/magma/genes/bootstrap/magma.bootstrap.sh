#!/bin/bash

if [[ ! $# -eq 2 ]]; then
  echo "Usage: magma.bootstrap.sh <data_dir> <destination_dir>"
  exit 1
fi
DATA_DIR=$1
INPUT_DIR=$2

aws s3 cp ${DATA_DIR}/bin/magma/magma_v1.07bb_static.zip ${INPUT_DIR}/inputs/
unzip ${INPUT_DIR}/inputs/magma_v1.07bb_static.zip -d ${INPUT_DIR}/inputs/magma
rm ${INPUT_DIR}/inputs/magma_v1.07bb_static.zip
