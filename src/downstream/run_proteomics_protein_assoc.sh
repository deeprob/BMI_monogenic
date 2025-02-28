#!/bin/bash

notebook_path="notebooks/bmi/14a_proteomics_protein_assoc.ipynb"


dx login --token TOKEN  

my_cmd="pip install openpyxl statsmodels scikit-learn && papermill 14a_proteomics_protein_assoc.ipynb 14a_proteomics_protein_assoc_out.ipynb"

dx run dxjupyterlab \
    -icmd="$my_cmd" \
    -iin="${notebook_path}" \
    -iduration=480 \
    --destination "notebooks/bmi/data/downstream/proteomics/" \
    --instance-type "mem3_ssd1_v2_x4" \
    --priority high -y

