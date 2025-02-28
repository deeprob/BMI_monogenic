#!/bin/bash

dx login --token TOKEN

# Define the range of values for c and lf
c_values=({1..22}) #(21 22)
ancestry=("eur" "sas" "afr" "eas" "amr" "mid" "oth") 

# Iterate over c and lf values
for c in "${c_values[@]}"; do
  echo $c
  for a in "${ancestry[@]}"; do
    echo $a
    # Call the helper script with c and lf as arguments
    ./2_run_regenie_step2_quant_helper.sh "$c" "$a"
  done
done
