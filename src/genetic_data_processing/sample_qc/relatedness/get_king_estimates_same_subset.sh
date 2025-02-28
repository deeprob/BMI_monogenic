#!/bin/bash

dx login --token TOKEN  

subset_num=$1
my_cmd="wget https://www.kingrelatedness.com/Linux-king.tar.gz && \
        tar -xzvf Linux-king.tar.gz && \
        ./king -b subset${subset_num}.bed --kinship --degree 3 --cpus 90 --prefix subset${subset_num}"

bed_file="/notebooks/wes/sample_qc/relatedness/subset${subset_num}.bed"
bim_file="/notebooks/wes/sample_qc/relatedness/subset${subset_num}.bim"
fam_file="/notebooks/wes/sample_qc/relatedness/subset${subset_num}.fam"


dx run app-swiss-army-knife \
    -icmd="${my_cmd}" \
    -iin="${bed_file}" \
    -iin="${bim_file}" \
    -iin="${fam_file}" \
    --destination "notebooks/wes/sample_qc/relatedness/" \
    --instance-type "mem2_ssd1_v2_x96" \
    --priority low -y

# Same subsets to run 1,2,3,4,5
# This script takes 22 mins to run and costs 0.26 pounds
