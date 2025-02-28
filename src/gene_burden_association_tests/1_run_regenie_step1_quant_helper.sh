#!/bin/bash

a=$1

bed_file="/notebooks/snp/liftover/ukb_c1-22_GRCh38_full_analysis_set_plus_decoy_hla_merged.bed"
bim_file="/notebooks/snp/liftover/ukb_c1-22_GRCh38_full_analysis_set_plus_decoy_hla_merged.bim"
fam_file="/notebooks/snp/liftover/ukb_c1-22_GRCh38_full_analysis_set_plus_decoy_hla_merged.fam"
snplist="/notebooks/snp/liftover_qc/ancestry/final_array_snps_GRCh38_qc_pass_pruned_${a}.snplist"
snpid="/notebooks/snp/liftover_qc/final_array_snps_GRCh38_qc_pass.id"
pheno_file="/notebooks/bmi/data/processed/${a}_phenotype.tsv.gz"

my_cmd="
    docker pull ghcr.io/rgcgithub/regenie/regenie:v3.5.gz
    docker run \
        --name regenie_run \
        -v "./:/proj_dir/" \
        ghcr.io/rgcgithub/regenie/regenie:v3.5.gz regenie \
        --step 1 \
        --bed /proj_dir/ukb_c1-22_GRCh38_full_analysis_set_plus_decoy_hla_merged \
        --extract /proj_dir/final_array_snps_GRCh38_qc_pass_pruned_${a}.snplist \
        --keep /proj_dir/final_array_snps_GRCh38_qc_pass.id \
        --phenoFile /proj_dir/${a}_phenotype.tsv.gz \
        --covarFile /proj_dir/${a}_phenotype.tsv.gz \
        --phenoColList bmi_rint,bmi \
        --covarColList age,age_2,age_sex,genetic_pca{1:10} \
        --catCovarList genetic_sex,exome_release_batch \
        --qt \
        --print-prs \
        --bsize 1000 \
        --lowmem \
        --lowmem-prefix /proj_dir/regenie_tmp_preds \
        --out /proj_dir/ukb_step1_quant \
    "

dx run app-swiss-army-knife \
    -icmd="$my_cmd" \
    -iin="${bed_file}" \
    -iin="${bim_file}" \
    -iin="${fam_file}" \
    -iin="${pheno_file}" \
    -iin="${snplist}" \
    -iin="${snpid}" \
    --destination "/notebooks/bmi/data/assoc/step1/quant/${a}/" \
    --instance-type "mem1_ssd1_v2_x36" \
    --priority low \
    -y
