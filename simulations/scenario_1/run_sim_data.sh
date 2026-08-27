#!/bin/bash

for dim in 4 28; do
  for replica in {1..6}; do
    echo "Running dimension $dim replica $replica"
    Rscript ./simulations/scenario_1/01__sim_data.R $dim $replica
  done
done
