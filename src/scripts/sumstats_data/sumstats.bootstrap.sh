#!/bin/bash

aws s3 cp s3://dig-ldsc-server/bin/snplist/ snplist/ --recursive

mkdir -p ./weights
for ANCESTRY in AFR AMR EAS EUR SAS
do
  mkdir -p ./weights/$ANCESTRY
  aws s3 cp s3://dig-ldsc-server/bin/weights/weights_$ANCESTRY.zip ./
  unzip weights_$ANCESTRY.zip -d ./weights/$ANCESTRY/
  rm weights_$ANCESTRY.zip
done
