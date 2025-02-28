dx login --token TOKEN

subset_num=$1
my_cmd="plink2 --bfile autosome_hqc --keep-col-match subsets.txt ${subset_num} \
        --make-bed --out subset${subset_num} --threads 32"

bed_file="/notebooks/wes/sample_qc/high_quality_variants/autosomes/autosome_hqc.bed"
bim_file="/notebooks/wes/sample_qc/high_quality_variants/autosomes/autosome_hqc.bim"
fam_file="/notebooks/wes/sample_qc/high_quality_variants/autosomes/autosome_hqc.fam"
subsets_file="/notebooks/wes/sample_qc/relatedness/subsets.txt"

dx run app-swiss-army-knife \
    -icmd="${my_cmd}" \
    -iin="${bed_file}" \
    -iin="${bim_file}" \
    -iin="${fam_file}" \
    -iin="${subsets_file}" \
    --destination "notebooks/wes/sample_qc/relatedness/" \
    --instance-type "mem2_ssd1_v2_x32" \
    --priority low -y
