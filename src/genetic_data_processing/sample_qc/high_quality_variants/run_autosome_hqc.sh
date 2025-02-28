#!/bin/bash

chr_num=$1

dx login --token TOKEN  

my_cmd="awk '{ \$6=\"NONE\"; print \$0 }' ukb23158_c${chr_num}_b0_v1.fam > temp.fam && mv temp.fam ukb23158_c${chr_num}_b0_v1.fam && \
        plink2 --bfile ukb23158_c${chr_num}_b0_v1 \
        --out chr${chr_num}_hqc \
        --maf 0.001 --hwe 1e-6 --geno 0.01 \
        --indep-pairwise 500 50 0.2 \
        --make-bed --write-snplist \
        --threads 32"

bed_file="/Bulk/Exome sequences/Population level exome OQFE variants, PLINK format - final release/ukb23158_c${chr_num}_b0_v1.bed"
bim_file="/Bulk/Exome sequences/Population level exome OQFE variants, PLINK format - final release/ukb23158_c${chr_num}_b0_v1.bim"
fam_file="/Bulk/Exome sequences/Population level exome OQFE variants, PLINK format - final release/ukb23158_c${chr_num}_b0_v1.fam"


dx run app-swiss-army-knife \
    -icmd="${my_cmd}" \
    -iin="${bed_file}" \
    -iin="${bim_file}" \
    -iin="${fam_file}" \
    --destination "notebooks/wes/sample_qc/high_quality_variants/chr${chr_num}" \
    --instance-type "mem2_ssd2_v2_x32" \
    --priority low -y
