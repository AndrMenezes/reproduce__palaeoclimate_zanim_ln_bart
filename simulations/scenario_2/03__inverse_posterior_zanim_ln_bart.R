rm(list = ls())
library(zanicc)
library(ggplot2)

args <- commandArgs(trailingOnly = TRUE)
d <- as.integer(args[1L])
p <- as.integer(args[2L])
replica <- as.integer(args[3L])

# d <- 28L
# p <- 3L
# replica <- 1L

# Path
path_local <- sprintf("./simulations/scenario_2/d=%i_p=%i", d, p)
path_data <- file.path(path_local, "data")
path_results <- file.path(path_local, "results", sprintf("replica=%i", replica))
forests_dir <- file.path(path_results, "forests")
if (!dir.exists(forests_dir)) dir.create(forests_dir, recursive = TRUE)

# Import data
list_data <- readRDS(file.path(path_data, sprintf("data_replica=%i.rds", replica)))
id_test <- list_data$id_test
n_test <- length(id_test)
Y_test <- list_data$Y[id_test, ]
X_test <- list_data$X[id_test, , drop = FALSE]

# Priors
mean_prior <- colMeans(list_data$X)
S_prior <- 1.5*cov(list_data$X)

# Load forward model
zanim_ln_bart <- load_model(model_dir = path_results)

# Generate uniform proposal
N_PROPOSAL <- 2000L
if (!file.exists(file.path(path_results, "x_proposal.rds"))) {
  set.seed(1212)
  x_proposal <- rconvexhull(n = N_PROPOSAL, X = list_data$X)
  saveRDS(x_proposal, file.path(path_results, "x_proposal.rds"))
} else {
  x_proposal <- readRDS(file.path(path_results, "x_proposal.rds"))
}

cat("SIR\n")

# SIR
ff_sir <- file.path(path_results, "sir.rds")
if (!file.exists(ff_sir)) {
  sir <- inverse_posterior_zanimlnbart(object = zanim_ln_bart, Y = Y_test,
                                       x_proposal = x_proposal,
                                       dir_posterior_fx = path_results,
                                       method = "sir")
  # Save results
  saveRDS(object = sir, file = ff_sir)
}
cat("ESS\n")

# Run eSS
ff_ess <- file.path(path_results, "ess.rds")
if (!file.exists(ff_ess)) {
  ess <- inverse_posterior_zanimlnbart(object = zanim_ln_bart, Y = Y_test,
                                       dir_posterior_fx = path_results,
                                       method = "ess",
                                       mean_prior = mean_prior,
                                       S_prior = S_prior,
                                       nburnin = 1L, n_particles = 200L)
  # Save results
  saveRDS(object = ess, file = ff_ess)
}

cat("cESS\n")
# H-representation for the cESS
hull <- geometry::convhulln(list_data$X, options = "Pp Fa Fx", output.options = "n")
normals <- hull$normals
A <- -normals[, -ncol(normals), drop = FALSE]
b <- -normals[, ncol(normals)]

dim(A)
length(b)

# Run ceSS
ff_cess <- file.path(path_results, "cess.rds")
if (!file.exists(ff_cess)) {
  cess <- inverse_posterior_zanimlnbart(object = zanim_ln_bart,
                                        Y = Y_test, method = "cess",
                                        mean_prior = mean_prior,
                                        S_prior = S_prior,
                                        Amat = A, bvec = b,
                                        eta = 500L,
                                        nburnin = 1L, n_particles = 200L)
  # Save results
  saveRDS(object = cess, file = ff_cess)
}
