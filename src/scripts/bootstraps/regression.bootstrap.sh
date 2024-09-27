#!/bin/bash

aws s3 cp s3://dig-ldsc-server/bin/regression_inputs/baseline_regression_inputs.zip baseline/
unzip baseline/baseline_regression_inputs.zip -d baseline/
rm baseline/baseline_regression_inputs.zip

# only EUR for now
aws s3 cp s3://dig-ldsc-server/bin/regression_inputs/tissue_regression_inputs_EUR.zip tissue/
unzip tissue/tissue_regression_inputs_EUR.zip -d tissue/
rm tissue/tissue_regression_inputs_EUR.zip

aws s cp s3://dig-ldsc-server/bin/regression_inputs/overlap.zip overlap/
unzip overlap/overlap.zip -d overlap/
rm overlap/overlap.zip
