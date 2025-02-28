#!/bin/bash

chr_num=$1
notebook_path="notebooks/wes/variant_annot/chr${chr_num}.ipynb"


dx login --token TOKEN  

my_cmd="papermill chr${chr_num}.ipynb chr${chr_num}_out.ipynb"

dx run dxjupyterlab_spark_cluster \
    -ifeature="HAIL-VEP" \
    -icmd="$my_cmd" \
    -iin="${notebook_path}" \
    -iduration=120 \
    --destination "notebooks/wes/variant_annot" \
    --instance-type "mem2_ssd1_v2_x32" \
    --instance-count 30 \
    --priority low -y
