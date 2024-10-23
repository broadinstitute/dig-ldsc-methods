#!/bin/bash

for ANCESTRY in AFR AMR EAS EUR SAS
do
  mkdir -p ./baseline/$ANCESTRY
  aws s3 cp s3://dig-ldsc-server/bin/baseline/baseline_$ANCESTRY.zip ./
  unzip baseline_$ANCESTRY.zip -d ./baseline/$ANCESTRY/
  rm baseline_$ANCESTRY.zip
done

for ANCESTRY in AFR AMR EAS EUR SAS
do
  mkdir -p ./tissue/$ANCESTRY
  aws s3 cp s3://dig-ldsc-server/bin/tissue/tissue_$ANCESTRY.zip ./
  unzip tissue_$ANCESTRY.zip -d ./tissue/$ANCESTRY/
  rm tissue_$ANCESTRY.zip
done

for ANCESTRY in AFR AMR EAS EUR SAS
do
  mkdir -p ./frq/$ANCESTRY
  aws s3 cp s3://dig-ldsc-server/bin/frq/frq_$ANCESTRY.zip ./
  unzip frq_$ANCESTRY.zip -d ./frq/$ANCESTRY/
  rm frq_$ANCESTRY.zip
done
