rm(list = ls())
library(zanicc)
source("./application/_cv_function.R")

# Settings
NTREES_THETA <- 100L
NTREES_ZETA <- 100L
NDPOST <- 4000L
NSKIP <- 5000L
NVALIDATION <- 1000L

# Paths
path_local <- "./application"
results_dir <- file.path(path_local, "results", "validation")
forests_dir <- file.path(results_dir, "forests")
if (!dir.exists(forests_dir)) dir.create(forests_dir, recursive = TRUE)

# Load data
data("pollen_climate")
Y <- pollen_climate$Y
X <- pollen_climate$X[, c("gdd5", "mtco", "aet.pet")]
X <- scale(X)
d <- ncol(Y)
n <- nrow(Y)

# Split data
if (!file.exists(file.path(results_dir, "idx_validation.rds"))) {
  set.seed(1212)
  idx_val <- sample.int(n = n, size = NVALIDATION)
  saveRDS(object = idx_val, file = file.path(results_dir, "idx_validation.rds"))
}
idx_val <- readRDS(file = file.path(results_dir, "idx_validation.rds"))
Y_val <- Y[idx_val, , drop = FALSE]
X_val <- X[idx_val, , drop = FALSE]
Y_train <- Y[-idx_val, , drop = FALSE]
X_train <- X[-idx_val, , drop = FALSE]


cat("Run forward model\n\n")
# Run forward model
if (!file.exists(file.path(results_dir, "mod.rds"))) {
  zanim_ln_bart <- zanicc(Y = Y_train, X_count = X_train, X_zi = X_train,
                          model = "zanim_ln_bart", ntrees_theta = NTREES_THETA,
                          ntrees_zeta = NTREES_ZETA, ndpost = NDPOST,
                          nskip = NSKIP, keep_draws = FALSE,
                          save_trees = TRUE, forests_dir = forests_dir)
  save_model(object = zanim_ln_bart, model_dir = results_dir)
}
zanim_ln_bart <- load_model(model_dir = results_dir)
zanim_ln_bart$elapsed_time[3]/60/60

# Proposal
N_PROPOSAL <- 4000L
if (!file.exists(file.path(results_dir, "x_proposal.rds"))) {
  x_proposal <- runifconvexhull(n = N_PROPOSAL, X = unique(X))
  saveRDS(x_proposal, file.path(results_dir, "x_proposal.rds"))
}
x_proposal <- readRDS(file.path(results_dir, "x_proposal.rds"))

# SIR
ff_sir <- file.path(results_dir, "sir.rds")
if (!file.exists(ff_sir)) {
  cat("\n\n SIR\n\n")
  sir <- inverse_posterior_zanimlnbart(object = zanim_ln_bart, Y = Y_val,
                                       x_proposal = x_proposal,
                                       dir_posterior_fx = forests_dir,
                                       method = "sir")
  saveRDS(object = sir, file = ff_sir)
}

# Priors for the ESS
S_prior <- 3.125*cov(X)
mean_prior <- colMeans(X)
# ESS
ff_ess <- file.path(results_dir, "ess.rds")
if (!file.exists(ff_ess)) {
  cat("\n\n ESS \n\n")
  ess <- inverse_posterior_zanimlnbart(object = zanim_ln_bart,
                                       Y = Y_val, method = "ess",
                                       mean_prior = mean_prior,
                                       S_prior = S_prior,
                                       nburnin = 1L,
                                       n_particles = 200L)
  saveRDS(object = ess, file = ff_ess)
}

# H-representation for the cESS
hull <- geometry::convhulln(X, options = "Pp Fa Fx", output.options = "n")
normals <- hull$normals
A <- -normals[, -ncol(normals), drop = FALSE]
b <- -normals[, ncol(normals)]
ff_cess <- file.path(results_dir, "cess.rds")
if (!file.exists(ff_cess)) {
  cat("\n\n cESS \n\n")
  cess <- inverse_posterior_zanimlnbart(object = zanim_ln_bart,
                                        Y = Y_val,
                                        method = "cess",
                                        mean_prior = mean_prior,
                                        S_prior = S_prior,
                                        Amat = A, bvec = b,
                                        eta = 500L,
                                        nburnin = 1L, n_particles = 200L)
  saveRDS(object = cess, file = ff_cess)
}


# Load posteriors
sir <- readRDS(ff_sir)
ess <- readRDS(ff_ess)
cess <- readRDS(ff_cess)

c(sir_sampling = attr(sir, "elapsed_time")[3],
  sir_predict = attr(sir, "elapsed_time_predict")[3],
  ess = attr(ess, "elapsed_time")[3],
  cess = attr(cess, "elapsed_time")[3])

format_sec_hour <- function(x) {
  sprintf("%dh %1dmin", x %/% 3600, (x %% 3600) %/% 60)
}

times <- c(sir = paste0(formatC(attr(sir, "elapsed_time")[3], digits = 1, format = "f"),
                        " (",
                        formatC(attr(sir, "elapsed_time_predict")[3], digits = 1, format = "f"),
                        ")"),
           ess = formatC(attr(ess, "elapsed_time")[3], digits = 1, format = "f") 
           ,
           cess = formatC(attr(cess, "elapsed_time")[3], digits = 1, format = "f"))
times <- c(
  sir = format_sec_hour(attr(sir, "elapsed_time")[3] + attr(sir, "elapsed_time_predict")[3]),
  ess = format_sec_hour(attr(ess, "elapsed_time")[3]),
  cess = format_sec_hour(attr(cess, "elapsed_time")[3]))

format_sec_hour(attr(sir, "elapsed_time_predict")[3])
format_sec_hour(attr(sir, "elapsed_time")[3])

# # Compute metrics
# (metrics_sir <- posterior_risk(x = X_val, draws = sir, joint_coverage = TRUE))
# (metrics_ess <- posterior_risk(x = X_val, draws = ess, joint_coverage = TRUE))
# (metrics_cess <- posterior_risk(x = X_val, draws = cess, joint_coverage = TRUE))
# 
# # Save
# saveRDS(object = metrics_sir, file = file.path(results_dir, "metrics_sir.rds"))
# saveRDS(object = metrics_ess, file = file.path(results_dir, "metrics_ess.rds"))
# saveRDS(object = metrics_cess, file = file.path(results_dir, "metrics_cess.rds"))

# Load
metrics_sir <- readRDS(file = file.path(results_dir, "metrics_sir.rds"))
metrics_ess <- readRDS(file = file.path(results_dir, "metrics_ess.rds"))
metrics_cess = readRDS(file = file.path(results_dir, "metrics_cess.rds"))

tab <- cbind(sir = metrics_sir, ess = metrics_ess, cess = metrics_cess)
chosen_rows <- c("es", "msep", "mae",  "coverage_95")
xtab <- data.frame(t(tab[chosen_rows, ])) 
xtab <- cbind(xtab, times = times)
rownames(xtab) <- c("SIR", "ESS", "cESS")
print(xtable::xtable(xtab, digits = 3, include.rownames = TRUE))
