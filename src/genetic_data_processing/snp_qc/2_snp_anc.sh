#!/bin/bash

dx login --token TOKEN

ancestry=("eur" "sas" "afr" "eas" "amr" "mid" "oth")

for a in "${ancestry[@]}"; do
    echo $a
    # Call the helper script with arguments
    ./0_snp_anc_helper.sh "$a"
done
