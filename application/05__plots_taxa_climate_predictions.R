rm(list = ls())
library(zanicc)
library(ggplot2)
library(cowplot)
theme_set(
  theme_cowplot(font_size = 16, font_family = "Palatino") +
    background_grid() +
    theme(legend.position = "top")
)

# Paths
path_local <- "./application"
data_dir <- file.path(path_local, "data")
results_dir <- file.path(path_local, "results", "forward_module")
forests_dir <- file.path(results_dir, "forests")
if (!dir.exists(forests_dir)) dir.create(forests_dir, recursive = TRUE)

# Load grid to perform the predictions
X_grid <- readRDS(file = file.path(data_dir, "X_3d_buffer.rds"))
head(X_grid)
n_grid <- nrow(X_grid)

# Load modern data
data("pollen_climate", package = "zanicc")
df_X <- data.frame(pollen_climate$X[, c("gdd5", "mtco", "aet.pet")])
pollen_names <- colnames(pollen_climate$Y)
d <- length(pollen_names)

# Load model
zanim_ln_bart <- load_model(model_dir = results_dir)

# Compute predictions and save them in disk
if (!file.exists(results_dir, "theta_ij.bin")) {
  X_grid_scaled <- scale(X_grid)
  predict(zanim_ln_bart, newdata = X_grid_scaled, type = "theta",
          output_dir = results_dir, load = FALSE)
  predict(zanim_ln_bart, newdata = X_grid_scaled, type = "zeta",
          output_dir = results_dir, load = FALSE)
}

# Load theta and zeta predictions
NDPOST <- 4000L
dtheta <- zanicc::load_bin_predictions(fname = file.path(results_dir, "theta_ij.bin"),
                                       n = n_grid, d = d, m = NDPOST)
dzeta <- zanicc::load_bin_predictions(fname = file.path(results_dir, "zeta_ij.bin"),
                                      n = n_grid, d = d, m = NDPOST)

# Compute and save the vartheta predictions
if (!file.exists(results_dir, "vartheta.rds")) {
  # Load draws of the chol(Sigma_V)
  dm1 <- d - 1L
  fname <- file.path(forests_dir, "chol_Sigma_V.bin")
  draws_chol_Sigma_V <- zanicc::load_bin_predictions(fname = fname, n = dm1,
                                                     d = dm1, m = NDPOST)
  Bt <- zanim_ln_bart$Bt

  NDPOST <- dim(dtheta)[3L]
  vartheta <- array(dim = dim(dtheta))
  # For each draw compute vartheta simulating from the model conditional on the
  # posterior draw of the latent variables
  for (k in seq_len(NDPOST)) {
    cat(k, "\n")
    # Generate the z's
    tmp <- lapply(seq_len(n_grid), function(i) {
      z <- stats::rbinom(n = d, size = 1, prob = 1.0 - dzeta[i,,k])
      is_zero <- z == 0L
      if (all(is_zero)) {
        vt <- rep(0.0, d)
      }
      else if (sum(is_zero) == d - 1L) {
        vt <- rep(0.0, d)
        vt[!is_zero] <- 1.0
      } else {
        v <- stats::rnorm(dm1) %*% draws_chol_Sigma_V[,,k]
        u <- drop(v %*% Bt)
        vt <- dtheta[i,,k] * z * exp(u)
        vt <- vt / sum(vt)
      }
      vt
    })
    vartheta[,,k] <- do.call(rbind, tmp)
  }
  # Saving
  saveRDS(object = vartheta, file = file.path(results_dir, "vartheta.rds"))
}
vartheta <- readRDS(file = file.path(results_dir, "vartheta.rds"))


# Computing the posterior summaries of the predictions -------------------------
data_summaries_theta <- zanicc::summarise_draws_3d(dtheta)
data_summaries_zeta <- zanicc::summarise_draws_3d(dzeta)
data_summaries_vartheta <- zanicc::summarise_draws_3d(vartheta)
data_summaries_theta$gdd5 <- data_summaries_zeta$gdd5 <- data_summaries_vartheta$gdd5 <- rep(X_grid[, "gdd5"], times = d)
data_summaries_theta$mtco <- data_summaries_zeta$mtco <- data_summaries_vartheta$mtco <- rep(X_grid[, "mtco"], times = d)
data_summaries_theta$aet.pet <- data_summaries_zeta$aet.pet <- data_summaries_vartheta$aet.pet <- rep(X_grid[, "aet.pet"], times = d)
data_summaries_theta$taxa <- data_summaries_zeta$taxa <- data_summaries_vartheta$taxa <- pollen_names[data_summaries_vartheta$category]


# Ex(im)port summaries ---------------------------------------------------------

# saveRDS(object = data_summaries_theta, file = file.path(results_dir, "summaries_3d_theta.rds"))
# saveRDS(object = data_summaries_zeta, file = file.path(results_dir, "summaries_3d_zeta.rds"))
# saveRDS(object = data_summaries_vartheta, file = file.path(results_dir, "summaries_3d_vartheta.rds"))

data_summaries_theta <- readRDS(file = file.path(results_dir, "summaries_3d_theta.rds"))
data_summaries_zeta <- readRDS(file = file.path(results_dir, "summaries_3d_zeta.rds"))
data_summaries_vartheta <- readRDS(file = file.path(results_dir, "summaries_3d_vartheta.rds"))


# Wrapper functions to plot ----------------------------------------------------
plot_bivariate_predictions <- function(data_summaries, df_X,
                                       x1 = "gdd5", x2 = "mtco",
                                       parameter = c("theta", "zeta", "vartheta"),
                                       breaks_ = NULL, nbreaks = 4L,
                                       full_caption = TRUE) {
  if (full_caption) {
    if (parameter == "zeta")
      string_cap <- r'(Population-level structural zero probabilities, $\zeta\!{}_{ij}$)'
    if (parameter == "theta")
      string_cap <- r'(Population-level compositional probabilities, $\theta\!{}_{ij}$)'
    if (parameter == "vartheta") {
      string_cap <- r'(Individual-level compositional probabilities, $\vartheta\!{}_{ij}$)'
    }
  } else {
    if (parameter == "zeta")
      string_cap <- r'($\zeta\!{}_{ij}$)'
    if (parameter == "theta")
      string_cap <- r'($\theta\!{}_{ij}$)'
    if (parameter == "vartheta") {
      string_cap <- r'($\vartheta\!{}_{ij}$)'
    }
  }
  if (parameter == "vartheta")
    data_summaries$median[data_summaries$median == 0] <- NA_real_

  label_x1 <- c("gdd5" = "GDD5 (°C days)", mtco = "MTCO (°C)",
                "aet.pet" = "AET/PET")[x1]
  label_x2 <- c("gdd5" = "GDD5 (°C days)", mtco = "MTCO (°C)",
                "aet.pet" = "AET/PET")[x2]

  if (is.null(breaks_)) {
    if (parameter == "zeta") {
      breaks_ <- unique(round(pretty(x = data_summaries$median, n = nbreaks), 3))
    } else {
      breaks_ <- pretty(x = log(data_summaries$median / (1.0 - data_summaries$median)),
                        n = nbreaks)
      breaks_ <- unique(round(1/(1 + exp(-breaks_)), 3))
    }
  }
  ggplot(data_summaries) +
    geom_point(data = df_X,
               aes(x = .data[[x1]], y = .data[[x2]]), alpha = 0.3, size = 0.5) +
    geom_raster(aes(x = .data[[x1]], y = .data[[x2]], fill = median), alpha = 0.8, interpolate = FALSE) +
    colorspace::scale_fill_continuous_sequential(palette = "Viridis",
                                                  transform = if (parameter == "zeta") "identity" else "logit",
                                                  rev = FALSE,
                                                  na.value = "black",
                                                  breaks = breaks_) +
    guides(fill = guide_colourbar(barwidth = 16, barheight = 0.75, title.position = "top")) +
    labs(x = label_x1, y = label_x2,
         fill = latex2exp::TeX(string_cap)) +
    scale_x_continuous(breaks = scales::pretty_breaks(8)) +
    scale_y_continuous(breaks = scales::pretty_breaks(8))
}

# Plot all bivariate pairs for all taxa ----------------------------------------
comb_pairs <- combn(x = c("gdd5", "mtco", "aet.pet"), m = 2)
list_plots <- list(gdd5_mtco = vector(mode = "list", length = d),
                   gdd5_aetpet = vector(mode = "list", length = d),
                   mtco_aetpet = vector(mode = "list", length = d))
for (j in seq_along(pollen_names)) {
  # j=16
  cat(j, "of", d, "\n")
  cur_taxa <- pollen_names[j]
  data_theta_cur <- dplyr::filter(data_summaries_theta, taxa == cur_taxa)
  data_zeta_cur <- dplyr::filter(data_summaries_zeta, taxa == cur_taxa)
  data_vartheta_cur <- dplyr::filter(data_summaries_vartheta, taxa == cur_taxa)
  for (k in seq_len(ncol(comb_pairs))) {
    full_cap <- if (k == 1) TRUE else FALSE
    v1 <- comb_pairs[1, k]
    v2 <- comb_pairs[2, k]
    p_theta <- plot_bivariate_predictions(data_summaries = data_theta_cur, df_X = df_X,
                                          x1 = v1, x2 = v2, parameter = "theta",
                                          full_caption = full_cap)
    p_zeta <- plot_bivariate_predictions(data_summaries = data_zeta_cur, df_X = df_X,
                                         x1 = v1, x2 = v2, parameter = "zeta",
                                         full_caption = full_cap)
    p_vartheta <- plot_bivariate_predictions(data_summaries = data_vartheta_cur,
                                             df_X = df_X, x1 = v1, x2 = v2,
                                             parameter = "vartheta",
                                             full_caption = full_cap)
    list_plots[[k]][[j]] <- plot_grid(p_theta, p_zeta, p_vartheta, nrow = 1)
  }
}

# Combine the three pairs of bivariate plots for all categories and save them
for (j in seq_along(pollen_names)) {
  cat(j, "of", d, "\n")
  cur_taxa <- pollen_names[j]
  ff <- paste0("climate_predictions__", tolower(cur_taxa), ".png")
  pp <- plot_grid(list_plots$gdd5_mtco[[j]], list_plots$gdd5_aetpet[[j]],
                  list_plots$mtco_aetpet[[j]], ncol = 1,
                  labels = "AUTO", align = "hv")
  save_plot(filename = file.path(results_dir, ff), plot = pp, base_height = 12)
}


# Plot the chosen categories  --------------------------------------------------

p_mtco_gdd5_chosen <- plot_grid(
  plot_bivariate_predictions(data_summaries = dplyr::filter(data_summaries_theta, taxa == "Juniperus"),
                             df_X = df_X, x1 = "gdd5", x2 = "mtco",
                             parameter = "theta", full_caption = FALSE),
  plot_bivariate_predictions(data_summaries = dplyr::filter(data_summaries_zeta, taxa == "Juniperus"),
                             df_X = df_X, x1 = "gdd5", x2 = "mtco",
                             parameter = "zeta", full_caption = FALSE),
  plot_bivariate_predictions(data_summaries = dplyr::filter(data_summaries_vartheta, taxa == "Juniperus"),
                             df_X = df_X, x1 = "gdd5", x2 = "mtco",
                             parameter = "vartheta", full_caption = FALSE),
  nrow = 1
)
title <- ggdraw() + draw_label(expression(italic(Juniperus)),
                               fontfamily = "Palatino", x = 0, hjust = 0) +
  theme(plot.margin = margin(0, 0, 0, 58))
p_mtco_gdd5_chosen <- plot_grid(title, p_mtco_gdd5_chosen, ncol = 1, rel_heights = c(0.1, 1))

p_gdd5_aetpet_chosen <- plot_grid(
  plot_bivariate_predictions(data_summaries = dplyr::filter(data_summaries_theta, taxa == "Pinus.D"),
                             df_X = df_X, x1 = "gdd5", x2 = "aet.pet",
                             parameter = "theta", full_caption = TRUE),
  plot_bivariate_predictions(data_summaries = dplyr::filter(data_summaries_zeta, taxa == "Pinus.D"),
                             df_X = df_X, x1 = "gdd5", x2 = "aet.pet",
                             parameter = "zeta", full_caption = TRUE),
  plot_bivariate_predictions(data_summaries = dplyr::filter(data_summaries_vartheta, taxa == "Pinus.D"),
                             df_X = df_X, x1 = "gdd5", x2 = "aet.pet",
                             parameter = "vartheta", full_caption = TRUE),
  nrow = 1
)
title <- ggdraw() + draw_label(expression(italic(Pinus)*.D), fontfamily = "Palatino", x = 0, hjust = 0) +
  theme(plot.margin = margin(0, 0, 0, 50))
p_gdd5_aetpet_chosen <- plot_grid(title, p_gdd5_aetpet_chosen, ncol = 1, rel_heights = c(0.1, 1))
p_gdd5_aetpet_chosen

p_mtco_aetpet_chosen <- plot_grid(
  plot_bivariate_predictions(data_summaries = dplyr::filter(data_summaries_theta, taxa == "Betula"),
                             df_X = df_X, x1 = "mtco", x2 = "aet.pet",
                             parameter = "theta", full_caption = FALSE),
  plot_bivariate_predictions(data_summaries = dplyr::filter(data_summaries_zeta, taxa == "Betula"),
                             df_X = df_X, x1 = "mtco", x2 = "aet.pet",
                             parameter = "zeta", full_caption = FALSE),
  plot_bivariate_predictions(data_summaries = dplyr::filter(data_summaries_vartheta, taxa == "Betula"),
                             df_X = df_X, x1 = "mtco", x2 = "aet.pet",
                             parameter = "vartheta", full_caption = FALSE),
  nrow = 1
)
title <- ggdraw() + draw_label(expression(italic(Betula)), fontfamily = "Palatino", x = 0, hjust = 0) +
  theme(plot.margin = margin(0, 0, 0, 50))
p_mtco_aetpet_chosen <- plot_grid(title, p_mtco_aetpet_chosen, ncol = 1, rel_heights = c(0.1, 1))

pp <- plot_grid(p_gdd5_aetpet_chosen, p_mtco_aetpet_chosen, p_mtco_gdd5_chosen,
                ncol = 1, labels = "AUTO", align = "hv")
save_plot(filename = file.path(results_dir, "climate_predictions_chosen_pinusd_betula_juniperus.png"),
          plot = pp, base_height = 12)
