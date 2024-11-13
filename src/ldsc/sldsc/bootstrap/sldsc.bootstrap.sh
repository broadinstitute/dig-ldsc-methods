#!/bin/bash

if [[ ! $# -eq 2 ]]; then
  echo "Usage: sldsc.bootstrap.sh <destination_dir> <ancestry>"
  exit 1
fi
DIR=$1
ANCESTRY=$2

aws s3 cp s3://dig-ldsc-server/bin/sldsc_inputs/sldsc_inputs.${ANCESTRY}.zip ${DIR}/inputs/
