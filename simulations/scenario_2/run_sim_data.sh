#!/bin/bash
for d in 28; do
  for p in 2 3; do
    for rep in {1..6}; do
      echo "Running d $d p $p replica $rep"
      Rscript ./simulations/scenario_2/01__sim_data.R $d $p $rep
    done
  done
done
