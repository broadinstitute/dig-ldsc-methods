#!/bin/bash

for ANCESTRY in AFR AMR EAS EUR SAS
do
  aws s3 cp s3://dig-ldsc-server/bin/sldsc_inputs/sldsc_inputs.${ANCESTRY}.zip inputs/
  unzip inputs/sldsc_inputs.${ANCESTRY}.zip -d inputs/${ANCESTRY}/
  rm inputs/sldsc_inputs.${ANCESTRY}.zip
done

for ANCESTRY in AFR AMR EAS EUR SAS
do
  mkdir -p ./weights/$ANCESTRY
  aws s3 cp s3://dig-ldsc-server/bin/weights/weights_$ANCESTRY.zip ./
  unzip weights_$ANCESTRY.zip -d ./weights/$ANCESTRY/
  rm weights_$ANCESTRY.zip
done
