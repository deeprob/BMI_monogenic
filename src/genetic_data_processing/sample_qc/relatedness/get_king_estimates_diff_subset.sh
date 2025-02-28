#!/bin/bash

dx login --token TOKEN  

subset_num1=$1
subset_num2=$2

my_cmd="wget https://www.kingrelatedness.com/Linux-king.tar.gz && \
        tar -xzvf Linux-king.tar.gz && \
        ./king -b subset${subset_num1}.bed,subset${subset_num2}.bed \
        --kinship --proj 93967 --degree 3 --cpus 90 --prefix subset${subset_num1}${subset_num2}"

bed_file1="/notebooks/wes/sample_qc/relatedness/subset${subset_num1}.bed"
bim_file1="/notebooks/wes/sample_qc/relatedness/subset${subset_num1}.bim"
fam_file1="/notebooks/wes/sample_qc/relatedness/subset${subset_num1}.fam"
bed_file2="/notebooks/wes/sample_qc/relatedness/subset${subset_num2}.bed"
bim_file2="/notebooks/wes/sample_qc/relatedness/subset${subset_num2}.bim"
fam_file2="/notebooks/wes/sample_qc/relatedness/subset${subset_num2}.fam"


dx run app-swiss-army-knife \
    -icmd="${my_cmd}" \
    -iin="${bed_file1}" \
    -iin="${bim_file1}" \
    -iin="${fam_file1}" \
    -iin="${bed_file2}" \
    -iin="${bim_file2}" \
    -iin="${fam_file2}" \
    --destination "notebooks/wes/sample_qc/relatedness/" \
    --instance-type "mem2_ssd1_v2_x96" \
    --priority low -y

# Different subsets to run 12,13,14,15,23,24,25,34,35,45
# This script takes 41 mins to run and costs 0.48 pounds
