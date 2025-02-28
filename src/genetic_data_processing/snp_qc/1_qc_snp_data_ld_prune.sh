#!/bin/bash
dx login --token TOKEN

my_cmd="plink2 --bfile ukb_c1-22_GRCh38_full_analysis_set_plus_decoy_hla_merged \
        --out final_array_snps_GRCh38_qc_pass_pruned \
        --extract final_array_snps_GRCh38_qc_pass.snplist \
        --indep-pairwise 1000 100 0.9 \
        --write-snplist \
        --threads 32"

bed_file="/notebooks/snp/liftover/ukb_c1-22_GRCh38_full_analysis_set_plus_decoy_hla_merged.bed"
bim_file="/notebooks/snp/liftover/ukb_c1-22_GRCh38_full_analysis_set_plus_decoy_hla_merged.bim"
fam_file="/notebooks/snp/liftover/ukb_c1-22_GRCh38_full_analysis_set_plus_decoy_hla_merged.fam"
vqc_file="/notebooks/snp/liftover_qc/final_array_snps_GRCh38_qc_pass.snplist"


dx run app-swiss-army-knife \
    -icmd="$my_cmd" \
    -iin="${bed_file}" \
    -iin="${bim_file}" \
    -iin="${fam_file}" \
    -iin="${vqc_file}" \
    --destination "/notebooks/snp/liftover_qc/" \
    --instance-type "mem2_ssd2_v2_x32" \
    -y
