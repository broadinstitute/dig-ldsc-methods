#!/bin/bash

if [[ ! $# -eq 3 ]]; then
  echo "Usage: g1000.bootstrap.sh <data_dir> <destination_dir> <ancestry>"
  exit 1
fi
DATA_DIR=$1
INPUT_DIR=$2
ANCESTRY=$3

mkdir -p ${INPUT_DIR}/g1000/${ANCESTRY}
aws s3 cp ${DATA_DIR}/bin/g1000/g1000_chr_${ANCESTRY}.zip ./
unzip g1000_chr_${ANCESTRY}.zip -d ${INPUT_DIR}/g1000/${ANCESTRY}
rm g1000_chr_${ANCESTRY}.zip
