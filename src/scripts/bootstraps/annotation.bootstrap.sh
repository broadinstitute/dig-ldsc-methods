#!/bin/bash

for ANCESTRY in AFR AMR EAS EUR SAS
do
  mkdir -p ./g1000/$ANCESTRY
  aws s3 cp s3://dig-ldsc-server/bin/g1000/g1000_chr_$ANCESTRY.zip ./
  unzip g1000_chr_$ANCESTRY.zip -d ./g1000/$ANCESTRY/
  rm g1000_chr_$ANCESTRY.zip
done
