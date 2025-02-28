#!/bin/bash

dx login --token TOKEN  

my_cmd="wget https://www.kingrelatedness.com/Linux-king.tar.gz && tar -xzvf Linux-king.tar.gz && \
        ./king -b autosome_hqc.bed --kinship --degree 2 --cpus 44"

bed_file="/notebooks/wes/sample_qc/high_quality_variants/autosomes/autosome_hqc.bed"
bim_file="/notebooks/wes/sample_qc/high_quality_variants/autosomes/autosome_hqc.bim"
fam_file="/notebooks/wes/sample_qc/high_quality_variants/autosomes/autosome_hqc.fam"


dx run app-swiss-army-knife \
    -icmd="${my_cmd}" \
    -iin="${bed_file}" \
    -iin="${bim_file}" \
    -iin="${fam_file}" \
    --destination "notebooks/wes/sample_qc/relatedness/" \
    --instance-type "mem1_ssd1_v2_x48" \
    --priority low -y

## This script takes 14 hours and costs 12.75 pounds
