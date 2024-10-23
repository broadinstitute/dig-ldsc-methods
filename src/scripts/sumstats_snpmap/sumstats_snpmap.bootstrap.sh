#!/bin/bash

for ancestry in AFR AMR EAS EUR SAS
do
  aws s3 cp s3://dig-ldsc-server/bin/g1000/g1000_chr_"$ancestry".zip data/g1000/"$ancestry"/
  unzip data/g1000/"$ancestry"/g1000_chr_"$ancestry".zip -d data/g1000/"$ancestry"/
  rm data/g1000/"$ancestry"/g1000_chr_"$ancestry".zip
done

aws s3 cp s3://dig-ldsc-server/bin/hapmap/w_hm3.snplist.bz2 .
bunzip2 w_hm3.snplist.bz2
mkdir -p data/hapmap
mv w_hm3.snplist data/hapmap/

aws s3 cp s3://dig-ldsc-server/bin/liftover/b37toHg38.over.chain liftover/
aws s3 cp s3://dig-ldsc-server/bin/liftover/liftOver.linux.x86_64.v287 liftover/
aws s3 cp s3://dig-ldsc-server/bin/liftover/liftOver.macOS.arm64 liftover/
