#!/bin/bash

aws s3 cp s3://dig-ldsc-server/bin/liftover/b37toHg38.over.chain liftover/
aws s3 cp s3://dig-ldsc-server/bin/liftover/liftOver.linux.x86_64.v287 liftover/
aws s3 cp s3://dig-ldsc-server/bin/liftover/liftOver.macOS.arm64 liftover/

aws s3 cp s3://dig-ldsc-server/bin/magma/dbSNP_common_GRCh37.csv files/
