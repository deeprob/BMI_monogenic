#!/bin/bash

chr_num="X"

dx login --token TOKEN  

my_cmd="plink2 --bfile chr${chr_num}_hqc \
        --chr X \
        --out chr${chr_num}_hqc_nopar --output-chr chrMT \
        --make-bed --threads 32"

bed_file="/notebooks/wes/sample_qc/high_quality_variants/chr${chr_num}/chr${chr_num}_hqc.bed"
bim_file="/notebooks/wes/sample_qc/high_quality_variants/chr${chr_num}/chr${chr_num}_hqc.bim"
fam_file="/notebooks/wes/sample_qc/high_quality_variants/chr${chr_num}/chr${chr_num}_hqc.fam"


dx run app-swiss-army-knife \
    -icmd="${my_cmd}" \
    -iin="${bed_file}" \
    -iin="${bim_file}" \
    -iin="${fam_file}" \
    --destination "notebooks/wes/sample_qc/high_quality_variants/chr${chr_num}" \
    --instance-type "mem2_ssd2_v2_x32" \
    --priority low -y
