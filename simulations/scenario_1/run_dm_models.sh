#!/bin/bash

for dim in 4 ; do
  for replica in {1..6}; do
    echo "Running dimension $dim replica $replica"
    Rscript ./simulations/scenario_1/04__dm_gp.R $dim $replica
    Rscript ./simulations/scenario_1/05__dm_bummer.R $dim $replica
  done
done
