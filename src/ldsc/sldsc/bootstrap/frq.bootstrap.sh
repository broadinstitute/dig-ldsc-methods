#!/bin/bash

if [[ ! $# -eq 3 ]]; then
  echo "Usage: frq.bootstrap.sh <data_dir> <destination_dir> <ancestry>"
  exit 1
fi
DATA_DIR=$1
INPUT_DIR=$2
ANCESTRY=$3

mkdir -p ${INPUT_DIR}/frq/${ANCESTRY}
aws s3 cp ${DATA_DIR}/bin/frq/frq_${ANCESTRY}.zip ./
unzip frq_${ANCESTRY}.zip -d ${INPUT_DIR}/frq/${ANCESTRY}
rm frq_${ANCESTRY}.zip
