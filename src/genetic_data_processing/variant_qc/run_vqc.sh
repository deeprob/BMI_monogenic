#!/bin/bash

chr_num=$1
notebook_path="notebooks/wes/variant_qc/0_initial_variant_qc_chr${chr_num}.ipynb"


dx login --token TOKEN  

my_cmd="pip install gnomad && papermill 0_initial_variant_qc_chr${chr_num}.ipynb 0_initial_variant_qc_chr${chr_num}_out.ipynb"

dx run dxjupyterlab_spark_cluster \
    -ifeature="HAIL" \
    -icmd="$my_cmd" \
    -iin="${notebook_path}" \
    -iduration=180 \
    --destination "notebooks/wes/variant_qc" \
    --instance-type "mem2_ssd1_v2_x32" \
    --instance-count 40 \
    --priority low -y
