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

if (!dir.exists(path_local)) dir.create(path_local, recursive = TRUE)
if (!dir.exists(path_data)) dir.create(path_data, recursive = TRUE)

# Function to simulate the data
sim_zanim_ln_gp <- function(n, d, n_trials, X_real, len_scale_theta = 2.0,
                            intercept_theta = stats::runif(d, -2.3, 2.3),
                            lo_pw = 0.01, up_pw = 0.2,
                            psi_range = c(0.30, 0.40)) {

  p <- ncol(X_real)
  if (p > 3L) stop("number of covariates should be 1, 2 or 3.")

  if (length(n_trials) != n) n_trials <- rep(n_trials[1L], n)

  # Simulate covariates
  if (p == 1L) X <- matrix(stats::runif(n, min(X_real), max(X_real)), ncol = 1L)
  else X <- zanicc::rconvexhull(n = n, X = unique(X_real))
  X <- sweep(X, 2, colMeans(X), "-")

  # Distance
  D <- as.matrix(fields::rdist(X))
  # Squared exponential kernel
  S_GP_theta <- exp(-D^2 / len_scale_theta) + diag(sqrt(.Machine$double.eps), nrow = n)
  chol_S_GP_theta <- chol(S_GP_theta)
  # S_GP_zeta <- exp(-D^2 / 10) + diag(sqrt(.Machine$double.eps), nrow = n)
  # chol_S_GP_zeta <- chol(S_GP_zeta)


  # Compositional
  fx_theta <- matrix(nrow = n, ncol = d)
  for (j in seq_len(d)) fx_theta[, j] <- drop(stats::rnorm(n) %*% chol_S_GP_theta)

  alphas <- exp(t(intercept_theta + t(fx_theta)))
  true_thetas <- sweep(alphas, 1, rowSums(alphas), "/")

  # Structural zero
  pw <- stats::runif(d, lo_pw, up_pw)#
  true_zetas <- matrix(nrow = n, ncol = d)
  for (j in seq_len(d)) true_zetas[, j] <- 1 - true_thetas[, j]^pw[j]

  # Random effect
  q_factors <- zanicc:::.ledermann(d)
  Gamma <- matrix(stats::runif(d * q_factors, 0, 1), d, q_factors)
  Psi <- diag(seq(psi_range[1L], psi_range[2L], length.out = d))
  Sigma_U <- tcrossprod(Gamma) + Psi
  chol_Sigma_U <- chol(Sigma_U)
  U <- matrix(0.0, nrow = n, ncol = d)
  for (i in seq_len(n)) U[i, ] <- drop(stats::rnorm(d) %*% chol_Sigma_U)


  # Simulate the counts
  Y <- Z <- matrix(0L, nrow = n, ncol = d)
  true_thetas <- true_varthetas <- matrix(0.0, nrow = n, ncol = d)
  z <- rep(1, d)
  for (i in seq_len(n)) {
    z <- stats::rbinom(d, 1L, prob = 1.0 - true_zetas[i, ])
    while (all(z == 0L)) z <- stats::rbinom(d, 1L, prob = 1 - true_zetas[i, ])
    Z[i, ] <- z
    is_zero <- z == 0L
    eU <- exp(U[i, ])
    true_varthetas[i, ] <- z * alphas[i, ] * eU / sum(z *  alphas[i, ] * eU)
    if (sum(is_zero) == d - 1L) {
      Y[i, ] <- rep(0L, d)
      Y[i, !is_zero] <- n_trials[i]
    } else {
      Y[i, ] <- stats::rmultinom(n = 1L, size = n_trials[i], prob = true_varthetas[i, ])
    }
  }
  print(cbind(structural_zeros = colMeans(1 - Z),
              sampling_zeros = colMeans(Y == 0) - colMeans(1 - Z)))

  data_sim <- data.frame(
    id = rep(seq_len(n), each = d),
    category = rep(seq_len(d), times = n),
    x1 = rep(X[, 1L], each = d),
    theta = c(t(true_thetas)),
    zeta = c(t(true_zetas)),
    vartheta = c(t(true_varthetas)),
    total = c(t(Y)),
    prop = c(apply(Y, 1L, function(y) y / sum(y))))
  if (p >= 2L) data_sim$x2 <- rep(X[, 2L], each = d)
  if (p == 3L) data_sim$x3 <- rep(X[, 3L], each = d)

  list_data <- list(df = data_sim, Y = Y, X = X, Z = Z, true_thetas = true_thetas,
                    true_zetas = true_zetas, true_varthetas = true_varthetas)
}

# Load pollen data
data("pollen_climate", package = "zanicc")
n_trials <- rowSums(pollen_climate$Y)
# prop.table(table(n_trials))

# Use the centered log-ratio of the empirical mean compositional as intercept
# for the compositional part
base_comp <- unname(colMeans(sweep(pollen_climate$Y, MARGIN = 1, n_trials, "/")))
intercept_theta <- log(base_comp / exp(mean(log(base_comp))))

# Seed
set.seed(replica)

# Sample size and number of categories
n_sample <- 1000L
n_trials <- sample(1000L:2000L, size = n_sample, replace = TRUE)

# Length scale of the GP associated to theta_ij
len_scale_theta <- 5.0

# Chose covariates
chosen <- switch(p, "mtco", c("gdd5", "mtco"), c("gdd5", "mtco", "aet.pet"))
X_real <- unique(scale(pollen_climate$X[, chosen, drop = FALSE]))

int_theta <- if (d == 28) intercept_theta else stats::runif(d, -2.3, 2.3)
list_data <- sim_zanim_ln_gp(n = n_sample, d = d, X_real = X_real,
                             n_trials = n_trials,
                             len_scale_theta = len_scale_theta,
                             intercept_theta = int_theta)
# Test sample
n_test <- 100L
list_data$id_test <- sample.int(n_sample, n_test)
# Save data
saveRDS(object = list_data,
        file = file.path(path_data, sprintf("data_replica=%i.rds", replica)))

# Plot of the true parameters
p_theta <- ggplot(data = list_data$df, aes(x = x1, y = x2)) +
  facet_wrap(~category) +
  geom_point(aes(col = theta), shape = 4) +
  scale_color_viridis_c(option = "C", limits = c(0, max(list_data$df$theta)))
p_zeta <- ggplot(data = list_data$df, aes(x = x1, y = x2)) +
  facet_wrap(~category) +
  geom_point(aes(col = zeta), shape = 4) +
  scale_color_viridis_c(option = "C", limits = c(0, max(list_data$df$zeta)))
p_vartheta <- ggplot(data = list_data$df, aes(x = x1, y = x2)) +
  facet_wrap(~category) +
  geom_point(aes(col = vartheta), shape = 4) +
  scale_color_viridis_c(option = "C", limits = c(0, max(list_data$df$vartheta)))
# Save plots
cowplot::save_plot(filename = file.path(path_data, sprintf("true_theta_replica=%i.png", replica)),
                   plot = p_theta, bg = "white", base_height = 8)
cowplot::save_plot(filename = file.path(path_data, sprintf("true_zeta_replica=%i.png", replica)),
                   plot = p_zeta, bg = "white", base_height = 8)
cowplot::save_plot(filename = file.path(path_data, sprintf("true_vartheta_replica=%i.png", replica)),
                   plot = p_vartheta, bg = "white", base_height = 8)
