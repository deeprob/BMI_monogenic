#!/bin/bash

dx login --token TOKEN  

my_cmd="for i in {1..22}; do echo \"chr\${i}_hqc\"; done > files_to_merge.txt && \
        plink2 --pmerge-list files_to_merge.txt bfile --make-bed \
        --out autosome_hqc \
        --threads 32"

bed_file1="/notebooks/wes/sample_qc/high_quality_variants/chr1/chr1_hqc.bed"
bim_file1="/notebooks/wes/sample_qc/high_quality_variants/chr1/chr1_hqc.bim"
fam_file1="/notebooks/wes/sample_qc/high_quality_variants/chr1/chr1_hqc.fam"
bed_file2="/notebooks/wes/sample_qc/high_quality_variants/chr2/chr2_hqc.bed"
bim_file2="/notebooks/wes/sample_qc/high_quality_variants/chr2/chr2_hqc.bim"
fam_file2="/notebooks/wes/sample_qc/high_quality_variants/chr2/chr2_hqc.fam"
bed_file3="/notebooks/wes/sample_qc/high_quality_variants/chr3/chr3_hqc.bed"
bim_file3="/notebooks/wes/sample_qc/high_quality_variants/chr3/chr3_hqc.bim"
fam_file3="/notebooks/wes/sample_qc/high_quality_variants/chr3/chr3_hqc.fam"
bed_file4="/notebooks/wes/sample_qc/high_quality_variants/chr4/chr4_hqc.bed"
bim_file4="/notebooks/wes/sample_qc/high_quality_variants/chr4/chr4_hqc.bim"
fam_file4="/notebooks/wes/sample_qc/high_quality_variants/chr4/chr4_hqc.fam"
bed_file5="/notebooks/wes/sample_qc/high_quality_variants/chr5/chr5_hqc.bed"
bim_file5="/notebooks/wes/sample_qc/high_quality_variants/chr5/chr5_hqc.bim"
fam_file5="/notebooks/wes/sample_qc/high_quality_variants/chr5/chr5_hqc.fam"
bed_file6="/notebooks/wes/sample_qc/high_quality_variants/chr6/chr6_hqc.bed"
bim_file6="/notebooks/wes/sample_qc/high_quality_variants/chr6/chr6_hqc.bim"
fam_file6="/notebooks/wes/sample_qc/high_quality_variants/chr6/chr6_hqc.fam"
bed_file7="/notebooks/wes/sample_qc/high_quality_variants/chr7/chr7_hqc.bed"
bim_file7="/notebooks/wes/sample_qc/high_quality_variants/chr7/chr7_hqc.bim"
fam_file7="/notebooks/wes/sample_qc/high_quality_variants/chr7/chr7_hqc.fam"
bed_file8="/notebooks/wes/sample_qc/high_quality_variants/chr8/chr8_hqc.bed"
bim_file8="/notebooks/wes/sample_qc/high_quality_variants/chr8/chr8_hqc.bim"
fam_file8="/notebooks/wes/sample_qc/high_quality_variants/chr8/chr8_hqc.fam"
bed_file9="/notebooks/wes/sample_qc/high_quality_variants/chr9/chr9_hqc.bed"
bim_file9="/notebooks/wes/sample_qc/high_quality_variants/chr9/chr9_hqc.bim"
fam_file9="/notebooks/wes/sample_qc/high_quality_variants/chr9/chr9_hqc.fam"
bed_file10="/notebooks/wes/sample_qc/high_quality_variants/chr10/chr10_hqc.bed"
bim_file10="/notebooks/wes/sample_qc/high_quality_variants/chr10/chr10_hqc.bim"
fam_file10="/notebooks/wes/sample_qc/high_quality_variants/chr10/chr10_hqc.fam"
bed_file11="/notebooks/wes/sample_qc/high_quality_variants/chr11/chr11_hqc.bed"
bim_file11="/notebooks/wes/sample_qc/high_quality_variants/chr11/chr11_hqc.bim"
fam_file11="/notebooks/wes/sample_qc/high_quality_variants/chr11/chr11_hqc.fam"
bed_file12="/notebooks/wes/sample_qc/high_quality_variants/chr12/chr12_hqc.bed"
bim_file12="/notebooks/wes/sample_qc/high_quality_variants/chr12/chr12_hqc.bim"
fam_file12="/notebooks/wes/sample_qc/high_quality_variants/chr12/chr12_hqc.fam"
bed_file13="/notebooks/wes/sample_qc/high_quality_variants/chr13/chr13_hqc.bed"
bim_file13="/notebooks/wes/sample_qc/high_quality_variants/chr13/chr13_hqc.bim"
fam_file13="/notebooks/wes/sample_qc/high_quality_variants/chr13/chr13_hqc.fam"
bed_file14="/notebooks/wes/sample_qc/high_quality_variants/chr14/chr14_hqc.bed"
bim_file14="/notebooks/wes/sample_qc/high_quality_variants/chr14/chr14_hqc.bim"
fam_file14="/notebooks/wes/sample_qc/high_quality_variants/chr14/chr14_hqc.fam"
bed_file15="/notebooks/wes/sample_qc/high_quality_variants/chr15/chr15_hqc.bed"
bim_file15="/notebooks/wes/sample_qc/high_quality_variants/chr15/chr15_hqc.bim"
fam_file15="/notebooks/wes/sample_qc/high_quality_variants/chr15/chr15_hqc.fam"
bed_file16="/notebooks/wes/sample_qc/high_quality_variants/chr16/chr16_hqc.bed"
bim_file16="/notebooks/wes/sample_qc/high_quality_variants/chr16/chr16_hqc.bim"
fam_file16="/notebooks/wes/sample_qc/high_quality_variants/chr16/chr16_hqc.fam"
bed_file17="/notebooks/wes/sample_qc/high_quality_variants/chr17/chr17_hqc.bed"
bim_file17="/notebooks/wes/sample_qc/high_quality_variants/chr17/chr17_hqc.bim"
fam_file17="/notebooks/wes/sample_qc/high_quality_variants/chr17/chr17_hqc.fam"
bed_file18="/notebooks/wes/sample_qc/high_quality_variants/chr18/chr18_hqc.bed"
bim_file18="/notebooks/wes/sample_qc/high_quality_variants/chr18/chr18_hqc.bim"
fam_file18="/notebooks/wes/sample_qc/high_quality_variants/chr18/chr18_hqc.fam"
bed_file19="/notebooks/wes/sample_qc/high_quality_variants/chr19/chr19_hqc.bed"
bim_file19="/notebooks/wes/sample_qc/high_quality_variants/chr19/chr19_hqc.bim"
fam_file19="/notebooks/wes/sample_qc/high_quality_variants/chr19/chr19_hqc.fam"
bed_file20="/notebooks/wes/sample_qc/high_quality_variants/chr20/chr20_hqc.bed"
bim_file20="/notebooks/wes/sample_qc/high_quality_variants/chr20/chr20_hqc.bim"
fam_file20="/notebooks/wes/sample_qc/high_quality_variants/chr20/chr20_hqc.fam"
bed_file21="/notebooks/wes/sample_qc/high_quality_variants/chr21/chr21_hqc.bed"
bim_file21="/notebooks/wes/sample_qc/high_quality_variants/chr21/chr21_hqc.bim"
fam_file21="/notebooks/wes/sample_qc/high_quality_variants/chr21/chr21_hqc.fam"
bed_file22="/notebooks/wes/sample_qc/high_quality_variants/chr22/chr22_hqc.bed"
bim_file22="/notebooks/wes/sample_qc/high_quality_variants/chr22/chr22_hqc.bim"
fam_file22="/notebooks/wes/sample_qc/high_quality_variants/chr22/chr22_hqc.fam"


dx run app-swiss-army-knife \
    -icmd="${my_cmd}" \
    -iin="${bed_file1}" -iin="${bim_file1}" -iin="${fam_file1}" \
    -iin="${bed_file2}" -iin="${bim_file2}" -iin="${fam_file2}" \
    -iin="${bed_file3}" -iin="${bim_file3}" -iin="${fam_file3}" \
    -iin="${bed_file4}" -iin="${bim_file4}" -iin="${fam_file4}" \
    -iin="${bed_file5}" -iin="${bim_file5}" -iin="${fam_file5}" \
    -iin="${bed_file6}" -iin="${bim_file6}" -iin="${fam_file6}" \
    -iin="${bed_file7}" -iin="${bim_file7}" -iin="${fam_file7}" \
    -iin="${bed_file8}" -iin="${bim_file8}" -iin="${fam_file8}" \
    -iin="${bed_file9}" -iin="${bim_file9}" -iin="${fam_file9}" \
    -iin="${bed_file10}" -iin="${bim_file10}" -iin="${fam_file10}" \
    -iin="${bed_file11}" -iin="${bim_file11}" -iin="${fam_file11}" \
    -iin="${bed_file12}" -iin="${bim_file12}" -iin="${fam_file12}" \
    -iin="${bed_file13}" -iin="${bim_file13}" -iin="${fam_file13}" \
    -iin="${bed_file14}" -iin="${bim_file14}" -iin="${fam_file14}" \
    -iin="${bed_file15}" -iin="${bim_file15}" -iin="${fam_file15}" \
    -iin="${bed_file16}" -iin="${bim_file16}" -iin="${fam_file16}" \
    -iin="${bed_file17}" -iin="${bim_file17}" -iin="${fam_file17}" \
    -iin="${bed_file18}" -iin="${bim_file18}" -iin="${fam_file18}" \
    -iin="${bed_file19}" -iin="${bim_file19}" -iin="${fam_file19}" \
    -iin="${bed_file20}" -iin="${bim_file20}" -iin="${fam_file20}" \
    -iin="${bed_file21}" -iin="${bim_file21}" -iin="${fam_file21}" \
    -iin="${bed_file22}" -iin="${bim_file22}" -iin="${fam_file22}" \
    --destination "notebooks/wes/sample_qc/high_quality_variants/autosomes/" \
    --instance-type "mem2_ssd2_v2_x32" \
    --priority normal -y
