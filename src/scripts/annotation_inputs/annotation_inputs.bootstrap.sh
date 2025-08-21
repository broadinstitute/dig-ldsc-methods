#!/bin/bash

mkdir -p ./heritability/
aws s3 cp s3://dig-ldsc-server/bin/sumstats/heritability.zip ./
unzip heritability.zip -d ./heritability
rm heritability.zip

for ANCESTRY in EUR
do
  mkdir -p ./sumstats/$ANCESTRY
  aws s3 cp s3://dig-ldsc-server/bin/sumstats/sumstats_$ANCESTRY.zip ./
  unzip sumstats_$ANCESTRY.zip -d ./sumstats/$ANCESTRY
  rm sumstats_$ANCESTRY.zip
done

for ANCESTRY in EUR
do
  mkdir -p ./baseline/$ANCESTRY
  aws s3 cp s3://dig-ldsc-server/bin/baseline/baseline_$ANCESTRY.zip ./
  unzip baseline_$ANCESTRY.zip -d ./baseline/$ANCESTRY
  rm baseline_$ANCESTRY.zip
done

for ANCESTRY in EUR
do
  mkdir -p ./weights/${ANCESTRY}
  aws s3 cp s3://dig-ldsc-server/bin/weights/weights_${ANCESTRY}.zip ./
  unzip weights_${ANCESTRY}.zip -d ./weights/${ANCESTRY}
  rm weights_${ANCESTRY}.zip
done

for ANCESTRY in EUR
do
  mkdir -p ./frq/${ANCESTRY}
  aws s3 cp s3://dig-ldsc-server/bin/frq/frq_${ANCESTRY}.zip ./
  unzip frq_${ANCESTRY}.zip -d ./frq/${ANCESTRY}
  rm frq_${ANCESTRY}.zip
done

cp ../../ldsc/sldsc/weights.py .
cp ../../ldsc/sldsc/sumstats.py .
cp ../../ldsc/sldsc/xtx_xty.py .
