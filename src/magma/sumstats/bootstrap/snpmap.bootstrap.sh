#!/bin/bash

if [[ ! $# -eq 4 ]]; then
  echo "Usage: snpmap.bootstrap.sh <data_dir> <destination_dir> <build_type> <genome_build>"
  exit 1
fi
DATA_DIR=$1
INPUT_DIR=$2
BUILD_TYPE=$3
GENOME_BUILD=$4

aws s3 cp ${DATA_DIR}/bin/magma/snpmap/sumstats.${BUILD_TYPE}.${GENOME_BUILD}.snpmap ${INPUT_DIR}/snpmap/
