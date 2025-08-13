#!/bin/bash

mkdir -p files
mkdir -p processed

aws s3 cp s3://dig-ldsc-server/bin/magma/NCBI37.3.gene.loc files/
aws s3 cp s3://dig-ldsc-server/bin/magma/dbSNP_common_GRCh37.csv files/

aws s3 cp ${DATA_DIR}/bin/magma/magma_v1.07bb_static.zip ${INPUT_DIR}/inputs/
unzip ${INPUT_DIR}/inputs/magma_v1.07bb_static.zip -d ${INPUT_DIR}/inputs/magma-linux
rm ${INPUT_DIR}/inputs/magma_v1.07bb_static.zip

aws s3 cp ${DATA_DIR}/bin/magma/magma_v1.07bb_mac.zip ${INPUT_DIR}/inputs/
unzip ${INPUT_DIR}/inputs/magma_v1.07bb_mac.zip -d ${INPUT_DIR}/inputs/magma-osx
rm ${INPUT_DIR}/inputs/magma_v1.07bb_mac.zip
