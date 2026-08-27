#!/bin/bash
for rep in {1..6}; do
  echo "Running d 28 p 3 replica $rep"
  Rscript ./simulations/scenario_2/02__zanim_ln_bart.R 28 3 $rep
  Rscript ./simulations/scenario_2/03__inverse_posterior_zanim_ln_bart.R 28 3 $rep
done
