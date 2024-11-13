#!/bin/bash

if [[ ! $# -eq 2 ]]; then
  echo "Usage: weights.bootstrap.sh <destination_dir> <ancestry>"
  exit 1
fi
DIR=$1
ANCESTRY=$2

mkdir -p ${DIR}/weights/${ANCESTRY}
aws s3 cp s3://dig-ldsc-server/bin/weights/weights_${ANCESTRY}.zip ./
unzip weights_${ANCESTRY}.zip -d ${DIR}/weights/${ANCESTRY}/
rm weights_${ANCESTRY}.zip
