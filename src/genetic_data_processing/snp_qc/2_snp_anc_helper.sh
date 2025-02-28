#!/bin/bash
a=$1

my_cmd="plink2 --bfile ukb_c1-22_GRCh38_full_analysis_set_plus_decoy_hla_merged \
        --extract final_array_snps_GRCh38_qc_pass_pruned.prune.in \
        --keep ${a}_samples.id \
        --out final_array_snps_GRCh38_qc_pass_pruned_${a} \
        --mac 5 \
        --write-snplist \
        --threads 32"

bed_file="/notebooks/snp/liftover/ukb_c1-22_GRCh38_full_analysis_set_plus_decoy_hla_merged.bed"
bim_file="/notebooks/snp/liftover/ukb_c1-22_GRCh38_full_analysis_set_plus_decoy_hla_merged.bim"
fam_file="/notebooks/snp/liftover/ukb_c1-22_GRCh38_full_analysis_set_plus_decoy_hla_merged.fam"
snplist="/notebooks/snp/liftover_qc/final_array_snps_GRCh38_qc_pass_pruned.prune.in"
sampleids="/notebooks/snp/liftover_qc/ancestry/${a}_samples.id"

dx run app-swiss-army-knife \
    -icmd="$my_cmd" \
    -iin="${bed_file}" \
    -iin="${bim_file}" \
    -iin="${fam_file}" \
    -iin="${snplist}" \
    -iin="${sampleids}" \
    --destination "/notebooks/snp/liftover_qc/ancestry/" \
    --instance-type "mem1_ssd1_v2_x36" \
    --priority low \
    -y
