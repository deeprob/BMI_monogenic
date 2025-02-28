#!/bin/bash
c=$1
a=$2

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
        --phenoColList bmi_rint,bmi \
        --covarColList age,age_2,age_sex,genetic_pca{1:10} \
        --catCovarList genetic_sex,exome_release_batch \
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
pheno_file="/notebooks/bmi/data/processed/${a}_phenotype.tsv.gz"
set_file="/notebooks/wes/burden_preparation/data/ukb_sets.tsv.gz"
anno_file="/notebooks/wes/burden_preparation/data/ukb_annotations.tsv.gz"
aaf_file="/notebooks/wes/burden_preparation/data/ukb_aafs.tsv.gz"
mask_file="/notebooks/wes/burden_preparation/data/ukb_masks.tsv.gz"
pred_list="/notebooks/bmi/data/assoc/step1/quant/${a}/ukb_step1_quant_pred.list"
prs_list="/notebooks/bmi/data/assoc/step1/quant/${a}/ukb_step1_quant_prs.list"
bmi_loco="/notebooks/bmi/data/assoc/step1/quant/${a}/ukb_step1_quant_1.loco"
bmi_prs="/notebooks/bmi/data/assoc/step1/quant/${a}/ukb_step1_quant_1.prs"
bmi_rint_loco="/notebooks/bmi/data/assoc/step1/quant/${a}/ukb_step1_quant_2.loco"
bmi_rint_prs="/notebooks/bmi/data/assoc/step1/quant/${a}/ukb_step1_quant_2.prs"


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
    -iin="${bmi_rint_loco}" \
    -iin="${bmi_rint_prs}" \
    --destination "/notebooks/bmi/data/assoc/step2/monogenic/${a}/chrm${c}/" \
    --instance-type "mem1_ssd1_v2_x36" \
    --priority low \
    -y --brief
