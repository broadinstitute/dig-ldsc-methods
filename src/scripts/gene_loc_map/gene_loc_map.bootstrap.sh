#!/bin/bash

mkdir -p files
mkdir -p processed

aws s3 cp s3://dig-ldsc-server/bin/gene_loc/NCBI37.3.plink.gene.loc files/
aws s3 cp s3://dig-ldsc-server/bin/gene_loc/gencode.gene.map files/
