#!/bin/bash
c=$1
m=$2
a=$3

my_cmd="
    docker pull ghcr.io/rgcgithub/regenie/regenie:v3.5.gz
    docker run \
        --name regenie_run${c} \
        -v "./:/proj_dir/" \
        ghcr.io/rgcgithub/regenie/regenie:v3.5.gz regenie \
        --step 2 \
        --pred /proj_dir/ukb_step1_quant_pred.list \
        --bgen /proj_dir/ukb23159_c${c}_b0_v1.bgen \
        --ref-first \
        --sample /proj_dir/ukb23159_c${c}_b0_v1.sample \
        --phenoFile /proj_dir/${a}_phenotype.tsv.gz \
        --covarFile /proj_dir/${a}_phenotype.tsv.gz \
        --phenoColList bmi,hba1c_df,hdl,ldl_sf \
        --covarColList age,genetic_pca{1:10} \
        --catCovarList genetic_sex \
        --set-list /proj_dir/ukb_sets.tsv.gz \
        --anno-file /proj_dir/ukb_annotations.tsv.gz \
        --aaf-file /proj_dir/ukb_aafs.tsv.gz \
        --mask-def /proj_dir/ukb_masks.tsv.gz \
        --minMAC 1 \
        --aaf-bins 0.001 \
        --bsize 200 \
        --out /proj_dir/bmi_quant \
    "

bgen_file="/Bulk/Exome sequences/Population level exome OQFE variants, BGEN format - final release//ukb23159_c${c}_b0_v1.bgen"
bgenidx_file="/Bulk/Exome sequences/Population level exome OQFE variants, BGEN format - final release//ukb23159_c${c}_b0_v1.bgen.bgi"
sample_file="/Bulk/Exome sequences/Population level exome OQFE variants, BGEN format - final release//ukb23159_c${c}_b0_v1.sample"
pheno_file="/notebooks/regenie/data/${a}_phenotype.tsv.gz"
set_file="/notebooks/regenie/data/ukb_sets.tsv.gz"
anno_file="/notebooks/regenie/data/ukb_annotations.tsv.gz"
aaf_file="/notebooks/regenie/data/ukb_aafs.tsv.gz"
mask_file="/notebooks/regenie/data/ukb_masks.tsv.gz"
pred_list="/notebooks/regenie/data/step1/quant/${a}/ukb_step1_quant_pred.list"
prs_list="/notebooks/regenie/data/step1/quant/${a}/ukb_step1_quant_prs.list"
bmi_loco="/notebooks/regenie/data/step1/quant/${a}/ukb_step1_quant_1.loco"
bmi_prs="/notebooks/regenie/data/step1/quant/${a}/ukb_step1_quant_1.prs"
hba1c_df_loco="/notebooks/regenie/data/step1/quant/${a}/ukb_step1_quant_2.loco"
hba1c_df_prs="/notebooks/regenie/data/step1/quant/${a}/ukb_step1_quant_2.prs" 
hdl_loco="/notebooks/regenie/data/step1/quant/${a}/ukb_step1_quant_3.loco"
hdl_prs="/notebooks/regenie/data/step1/quant/${a}/ukb_step1_quant_3.prs"
ldl_sf_loco="/notebooks/regenie/data/step1/quant/${a}/ukb_step1_quant_4.loco"
ldl_sf_prs="/notebooks/regenie/data/step1/quant/${a}/ukb_step1_quant_4.prs"


dx run app-swiss-army-knife \
    -icmd="$my_cmd" \
    -iin="${bgen_file}" \
    -iin="${bgenidx_file}" \
    -iin="${sample_file}" \
    -iin="${pheno_file}" \
    -iin="${set_file}" \
    -iin="${anno_file}" \
    -iin="${aaf_file}" \
    -iin="${mask_file}" \
    -iin="${pred_list}" \
    -iin="${prs_list}" \
    -iin="${bmi_loco}" \
    -iin="${bmi_prs}" \
    -iin="${hba1c_df_loco}" \
    -iin="${hba1c_df_prs}" \
    -iin="${hdl_loco}" \
    -iin="${hdl_prs}" \
    -iin="${ldl_sf_loco}" \
    -iin="${ldl_sf_prs}" \
    --destination "/notebooks/regenie/data/step2/monogenic/${a}/chrm${c}/output/" \
    --instance-type "mem2_ssd2_v2_x32" \
    -y --brief
