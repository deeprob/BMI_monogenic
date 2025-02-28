#!/bin/bash
dx login --token TOKEN

my_cmd="plink2 --bfile ukb_c1-22_GRCh38_full_analysis_set_plus_decoy_hla_merged \
        --out final_array_snps_GRCh38_qc_pass \
        --mac 100 --maf 0.01 --hwe 1e-15 --mind 0.1 --geno 0.1 \
        --write-snplist --write-samples --no-id-header \
        --threads 32"

bed_file="/notebooks/snp/liftover/ukb_c1-22_GRCh38_full_analysis_set_plus_decoy_hla_merged.bed"
bim_file="/notebooks/snp/liftover/ukb_c1-22_GRCh38_full_analysis_set_plus_decoy_hla_merged.bim"
fam_file="/notebooks/snp/liftover/ukb_c1-22_GRCh38_full_analysis_set_plus_decoy_hla_merged.fam"


dx run app-swiss-army-knife \
    -icmd="$my_cmd" \
    -iin="${bed_file}" \
    -iin="${bim_file}" \
    -iin="${fam_file}" \
    --destination "/notebooks/snp/liftover_qc/" \
    --instance-type "mem2_ssd2_v2_x32" \
    -y
