#!/bin/bash

if [[ ! $# -eq 4 ]]; then
  echo "Usage: snpmap.bootstrap.sh <destination_dir> <build_type> <genome_build> <ancestry>"
  exit 1
fi
DIR=$1
BUILD_TYPE=$2
GENOME_BUILD=$3
ANCESTRY=$4

aws s3 cp s3://dig-ldsc-server/bin/snpmap/sumstats.${BUILD_TYPE}.${GENOME_BUILD}.${ANCESTRY}.snpmap ${DIR}/snpmap/
