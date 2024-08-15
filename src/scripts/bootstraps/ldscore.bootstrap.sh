#!/bin/bash

aws s3 cp s3://dig-ldsc-server/bin/hapmap/hapmap3_snps.tgz .
tar -xf hapmap3_snps.tgz
mkdir -p hapmap
mv hapmap3_snps/* hapmap/
rm -r hapmap3_snps
rm hapmap3_snps.tgz

for ANCESTRY in AFR AMR EAS EUR SAS
do
  mkdir -p ./g1000/$ANCESTRY
  aws s3 cp s3://dig-ldsc-server/bin/g1000/g1000_chr_$ANCESTRY.zip ./
  unzip g1000_chr_$ANCESTRY.zip -d ./g1000/$ANCESTRY/
  rm g1000_chr_$ANCESTRY.zip
done
