#!/bin/bash

dx login --token TOKEN

ancestry=("eur" "sas" "afr" "eas" "amr" "mid" "oth")

for a in "${ancestry[@]}"; do
    echo $a
    # Call the helper script with arguments
    ./1_run_regenie_step1_quant_helper.sh "$a"
done
