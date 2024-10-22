#!/bin/bash

dx login --token TOKEN

# Define the range of values for c and lf
c_values=({1..22}) #(21 22)
mask=("PTV_Missense_strict") #"PTV" "PTV_Missense_lenient"
ancestry=("nonbritish") #"british" 

# Iterate over c and lf values
for m in "${mask[@]}"; do
  for c in "${c_values[@]}"; do
    echo $c
    for a in "${ancestry[@]}"; do
      echo $a
      # Call the helper script with c and lf as arguments
      # ./run_regenie_step2_binary_helper.sh "$c" "$m"
      ./run_regenie_step2_quant_helper.sh "$c" "$m" "$a"
    done
  done
done
