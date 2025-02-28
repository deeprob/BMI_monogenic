#!/bin/bash

chr_num="X"

dx login --token TOKEN  

my_cmd="plink2 --bfile chr${chr_num}_hqc_nopar \
        --out chr${chr_num}_hqc_nopar_pruned --output-chr chrMT \
        --exclude chr${chr_num}_hqc.prune.out \
        --make-bed \
        --threads 16"

bed_file="/notebooks/wes/sample_qc/high_quality_variants/chr${chr_num}/chr${chr_num}_hqc_nopar.bed"
bim_file="/notebooks/wes/sample_qc/high_quality_variants/chr${chr_num}/chr${chr_num}_hqc_nopar.bim"
fam_file="/notebooks/wes/sample_qc/high_quality_variants/chr${chr_num}/chr${chr_num}_hqc_nopar.fam"
exclude_file="/notebooks/wes/sample_qc/high_quality_variants/chr${chr_num}/chr${chr_num}_hqc.prune.out"


dx run app-swiss-army-knife \
    -icmd="${my_cmd}" \
    -iin="${bed_file}" \
    -iin="${bim_file}" \
    -iin="${fam_file}" \
    -iin="${exclude_file}" \
    --destination "notebooks/wes/sample_qc/high_quality_variants/chr${chr_num}" \
    --instance-type "mem1_ssd1_v2_x16" \
    --priority low -y
