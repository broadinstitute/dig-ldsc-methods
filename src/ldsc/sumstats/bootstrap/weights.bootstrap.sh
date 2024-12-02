#!/bin/bash

if [[ ! $# -eq 3 ]]; then
  echo "Usage: weights.bootstrap.sh <data_dir> <destination_dir> <ancestry>"
  exit 1
fi
DATA_DIR=$1
INPUT_DIR=$2
ANCESTRY=$3

mkdir -p ${INPUT_DIR}/weights/${ANCESTRY}
aws s3 cp ${DATA_DIR}/bin/weights/weights_${ANCESTRY}.zip ./
unzip weights_${ANCESTRY}.zip -d ${INPUT_DIR}/weights/${ANCESTRY}/
rm weights_${ANCESTRY}.zip
