rm(list = ls())
library(ggplot2)

args <- commandArgs(trailingOnly = TRUE)
d <- as.integer(args[1L])
replica <- as.integer(args[2L])

# d <- 4L
# replica <- 1L

path_local <- sprintf("./simulations/scenario_1/d=%i", d)
path_data <- file.path(path_local, "data")
path_results <- file.path(path_local, "results", sprintf("replica=%i", replica), "dm_gp")
if (!dir.exists(path_results)) dir.create(path_results, recursive = TRUE)

# Import data
list_data <- readRDS(file.path(path_data, sprintf("data_replica=%i.rds", replica)))
id_test <- list_data$id_test
n_test <- length(id_test)
Y_test <- list_data$Y[id_test, ]
X_test <- list_data$X[id_test, , drop = FALSE]
Y_train <- list_data$Y[-id_test, ]
X_train <- list_data$X[-id_test, , drop = FALSE]
true_thetas <- list_data$true_thetas[-id_test, ]
true_varthetas <- list_data$true_varthetas[-id_test, ]
true_zetas <- list_data$true_zetas[-id_test, ]
d <- ncol(Y_train)

# Define the model parameters
n_knots <- 30
X_knots <- seq(min(X_train, na.rm = TRUE) - 1.25*sd(X_train, na.rm = TRUE),
               max(X_train, na.rm = TRUE) + 1.25*sd(X_train, na.rm = TRUE),
               length = n_knots)
params <- list(
  n_adapt                    = 40000,
  n_mcmc                     = 10000,
  n_thin                     = 10,
  correlation_function       = "exponential",
  likelihood                 = "dirichlet-multinomial",
  function_type              = "gaussian-process",
  multiplicative_correlation = TRUE,
  ## prior on covariance to improve estimation
  eta                        = 8,
  additive_correlation       = FALSE,
  n_chains                   = 2,
  n_cores                    = 2,
  n_knots                    = n_knots,
  X_knots                    = X_knots)

# Fit
if (!file.exists(file.path(path_results, "mod.rds"))) {
  dm_gp <- BayesComposition::fit_compositional_data(y = Y_train, X = X_train, params = params,
                                                    progress_directory = path_results,
                                                    progress_file = "progress_file.txt")
  saveRDS(object = dm_gp, file = file.path(path_results, "mod.rds"))
}
dm_gp <- readRDS(file = file.path(path_results, "mod.rds"))

# Extract samples
samples <- BayesComposition::extract_compositional_samples(dm_gp)
alpha_post <- samples$alpha
theta_post <- sweep(alpha_post, MARGIN = c(1, 2),
                    STATS = apply(alpha_post, MARGIN = c(1, 2), sum),
                    FUN = "/")
theta_post <- aperm(theta_post, perm = c(2, 3, 1))
dim(theta_post)

# pdf(file.path(path_results, "diagnostics.pdf"), width = 8, height = 6)
par(mfrow = c(2, 2), mar = c(3, 3, 1, 1))
for (j in seq_len(d)) {
  plot(true_thetas[,j], rowMeans(theta_post[,j,]),
       xlab = "true", ylab = "estimate", main = sprintf("theta_{i%d}", j))
  abline(0, 1)
}

# Plot parameter against covariates
data_sim <- list_data$df[!(list_data$df$id %in% id_test), ]
data_sim$id <- rep(seq_len(nrow(Y_train)), each = d)
data_theta <- zanicc::summarise_draws_3d(x = theta_post)
data_theta$x <- rep(c(X_train), times = d)

saveRDS(object = data_theta, file = file.path(path_results, "posterior_theta.rds"))

p_theta <- ggplot(data = data_sim) +
  geom_line(mapping = aes(x = x, y = theta, col = "Truth", fill = "Truth"),
            linewidth = 0.8) +
  facet_wrap(~category, scales = "free_y") +
  geom_rug(data = dplyr::filter(data_sim, total == 0L),
           mapping = aes(y = NA_real_, x = x)) +
  geom_line(data = data_theta, mapping = aes(x = x, y = median),
            col = "dodgerblue") +
  geom_ribbon(data = data_theta,
              aes(x = x, ymin = ci_lower, ymax = ci_upper), fill = "dodgerblue",
              alpha = 0.3)
cowplot::save_plot(filename = file.path(path_results, "posterior_theta.png"),
                   plot = p_theta, bg = "white", base_height = 9)

# Posterior-predictive checks
n_trials <- rowSums(Y_train)
mc <- dim(alpha_post)[1L]
y_ppc <- array(dim = dim(alpha_post))
for (k in seq_len(mc)) {
  if (k %% 10 == 0L) cat(k, "\n")
  for (i in seq_len(nrow(Y_train))) {
    ld <- stats::rgamma(n = d, shape = alpha_post[k, i, ])
    y_ppc[k, i, ] <- drop(stats::rmultinom(n = 1L, size = n_trials[i],
                                           prob = ld / sum(ld)))
  }
}
png(filename = file.path(path_results, "ppc.png"), units = "in", width = 10,
    height = 7, res = 300)
out_ppc <- zanicc::plot_ppc(Y = Y_train, Y_ppc = y_ppc, output = TRUE)
graphics.off()
png(filename = file.path(path_results, "qqplots.png"), units = "in", width = 10,
    height = 7, res = 300)
zanicc::plot_qqplots(Y = Y_train, Y_ppc = y_ppc, relative = TRUE)
graphics.off()

# Inverse posterior
ini <- proc.time()

X_ini <- matrix(nrow = n_test, ncol = 1)
for (i in seq_len(n_test)) X_ini[i, ] <- stats::rnorm(1, mean = mean(X_train), sd = sd(X_train))

res <- BayesComposition::predict_compositional_data(
  y_reconstruct = Y_test,
  X_calibrate = X_ini,
  params = params,
  samples = samples,
  progress_directory = paste0(path_results, "/"),
  progress_file = "progress__inverse_posterior.txt")
end <- proc.time() - ini

ip_dm_gp <- array(dim = c(nrow(res$X), 1L, ncol(res$X)))
for (i in seq_len(ncol(res$X))) ip_dm_gp[,1L,i] <- res$X[, i]

saveRDS(ip_dm_gp,
        file = file.path(path_results, "ip_dm_gp.rds"))
