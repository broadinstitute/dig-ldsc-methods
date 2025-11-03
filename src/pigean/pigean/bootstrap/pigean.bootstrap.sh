#!/bin/bash

if [[ ! $# -eq 2 ]]; then
  echo "Usage: pigean.bootstrap.sh <data_dir> <destination_dir>"
  exit 1
fi
DATA_DIR=$1
INPUT_DIR=$2

aws s3 cp s3://dig-ldsc-server/bin/pigean/ . --recursive
