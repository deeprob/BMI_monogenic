#!/bin/bash

chrnum=$1
notebook_path="notebooks/wes/sample_qc/repartition/chr${chrnum}.ipynb"


dx login --token TOKEN  

my_cmd="pip install gnomad && papermill chr${chrnum}.ipynb chr${chrnum}_out.ipynb"

dx run dxjupyterlab_spark_cluster \
    -ifeature="HAIL" \
    -icmd="$my_cmd" \
    -iin="${notebook_path}" \
    -iduration=480 \
    --destination "notebooks/wes/sample_qc/repartition/" \
    --instance-type "mem2_ssd1_v2_x32" \
    --instance-count 40 \
    --priority low -y

# To track spark instance: https://job-xxxx.dnanexus.cloud:8081/jobs/
# chr21: mem2_ssd1_v2_x16 - nodes
