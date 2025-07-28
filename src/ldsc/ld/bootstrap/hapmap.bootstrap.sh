#!/bin/bash

if [[ ! $# -eq 2 ]]; then
  echo "Usage: hapmap.bootstrap.sh <data_dir> <destination_dir>"
  exit 1
fi
DATA_DIR=$1
INPUT_DIR=$2

mkdir -p ${INPUT_DIR}/hapmap
aws s3 cp ${DATA_DIR}/bin/hapmap/hapmap3_snps.tgz ./
tar -xf hapmap3_snps.tgz
mv hapmap3_snps/* ${INPUT_DIR}/hapmap
rm -r hapmap3_snps
rm hapmap3_snps.tgz
