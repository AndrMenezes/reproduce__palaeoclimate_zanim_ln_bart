# Reproducible Code for "Bayesian palaeoclimate reconstruction from zero-inflated count-compositional pollen data: A case study of Lago Grande di Monticchio in southern Italy"

This repository is associated with the article
["Bayesian palaeoclimate reconstruction from zero-inflated count-compositional pollen data: A case study of Lago Grande di Monticchio in southern Italy"](), and provides reproducible code for all simulations and the
palaeoclimate reconstruction of Lago Grande di Monticchio presented in the paper.

The results depends on the R package [`zanicc`](https://github.com/AndrMenezes/zanicc),
which implements the Bayesian modular framework for pollen-based palaeoclimate
reconstruction.


To install the development version of `zanicc`, you can use:
``` r
remotes::install_github("AndrMenezes/zanicc")
```

The repository is organised as follows:

- `simulations`: Code for the simulation studies presented in Section 3, comparing the proposed framework with existing methods.
  - `scenario_1`: Comparison of proposed framework against the DM-BUMMER and DM-GP models in reconstructing $p=1$ past climate.
  - `scenario_2`: Comparison of different sampling schemes for reconstructing $p=2$ and $p=3$ past climates.
- `application`: Code for reproduce the applications presented in Section 4. The additional plots of the estimated pollen-climate relationships are given in [here](application/results/forward_module).
