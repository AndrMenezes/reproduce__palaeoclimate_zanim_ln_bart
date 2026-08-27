#!/bin/bash

for dim in 4; do
  for replica in {1..6}; do
    echo "Running dimension $dim replica $replica"
    Rscript ./simulations/scenario_1/02__zanim_ln_bart.R $dim $replica
    Rscript ./simulations/scenario_1/03__inverse_posterior_zanim_ln_bart.R $dim $replica
  done
done
