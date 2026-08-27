rm(list = ls())
library(zanicc)

# Paths
path_local <- "./application"
results_dir <- file.path(path_local, "results")
forward_module_dir <- file.path(path_local, "results", "forward_module")
monticchio_dir <- file.path(path_local, "results", "monticchio")
if (!dir.exists(monticchio_dir)) dir.create(monticchio_dir, recursive = TRUE)

# Load modern data
data("pollen_climate")
X <- pollen_climate$X[, c("gdd5", "mtco", "aet.pet")]
X <- scale(X)
d <- ncol(pollen_climate$Y)

# Load fossil counts
data("pollen_monticchio")
Y_fossil <- pollen_monticchio

# Load ZANIM-LN-BART model fitted on the modern data
zanim_ln_bart <- load_model(model_dir = forward_module_dir)

# Generate proposal
N_PROPOSAL <- 2000L
if (!file.exists(file.path(monticchio_dir, "x_proposal.rds"))) {
  set.seed(1212)
  x_proposal <- runifconvexhull(n = N_PROPOSAL, X = unique(X))
  saveRDS(x_proposal, file.path(monticchio_dir, "x_proposal.rds"))
}
x_proposal <- readRDS(file.path(monticchio_dir, "x_proposal.rds"))

# SIR
ff_sir <- file.path(monticchio_dir, "sir.rds")
if (!file.exists(ff_sir)) {
  cat("\n\nSIR\n\n")
  sir <- inverse_posterior_zanimlnbart(object = zanim_ln_bart, Y = Y_fossil,
                                       x_proposal = x_proposal,
                                       dir_posterior_fx = monticchio_dir,
                                       method = "sir")
  saveRDS(object = sir, file = ff_sir)
}
