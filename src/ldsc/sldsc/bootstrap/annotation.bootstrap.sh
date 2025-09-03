#!/bin/bash

if [[ ! $# -eq 3 ]]; then
  echo "Usage: annotation.bootstrap.sh <data_dir> <destination_dir> <ancestry>"
  exit 1
fi
DATA_DIR=$1
INPUT_DIR=$2
ANCESTRY=$3

aws s3 cp ${DATA_DIR}/bin/annotation_inputs/annotation_inputs.${ANCESTRY}.zip ${INPUT_DIR}/inputs/
