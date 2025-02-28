#!/bin/bash

notebook_path="notebooks/wes/sample_qc/sample_qc.ipynb"


dx login --token TOKEN  

my_cmd="pip install gnomad && papermill sample_qc.ipynb sample_qc_out.ipynb"

dx run dxjupyterlab_spark_cluster \
    -ifeature="HAIL" \
    -icmd="$my_cmd" \
    -iin="${notebook_path}" \
    -iduration=600 \
    --destination "notebooks/wes/sample_qc/data/" \
    --instance-type "mem2_ssd2_v2_x32" \
    --instance-count 40 \
    --priority low -y

# To track spark instance: https://job-xxxx.dnanexus.cloud:8081/jobs/
# 
