rm(list = ls())
library(zanicc)
library(ggplot2)
library(cowplot)
source("./application/plot_bivariate.R")
theme_set(
  theme_cowplot(font_size = 16, font_family = "Palatino") +
    background_grid()
)

# Paths
path_local <- "./application"
results_dir <- file.path(path_local, "results", "monticchio")
forests_dir <- file.path(results_dir, "forests")
if (!dir.exists(forests_dir)) dir.create(forests_dir, recursive = TRUE)

# Load modern data
data("pollen_climate")
chosen_covariates <- c("gdd5", "mtco", "aet.pet")
X_modern <- pollen_climate$X[, chosen_covariates]
X <- scale(X_modern)
p <- ncol(X_modern)

# Load fossil counts
data("pollen_monticchio")
Y_fossil <- pollen_monticchio
age <- attr(Y_fossil, "age")
age <- age / 1000L

# Load the proposal
set.seed(1212)
N_PROPOSAL <- 20000L
x_unif <- runifconvexhull(n = N_PROPOSAL, X = unique(X))

# Load posterior
list.files(results_dir)
sir <- readRDS(file = file.path(results_dir, "sir.rds"))
attr(sir, "elapsed_time")[3]
attr(sir, "elapsed_time_predict")[3]

# Re-scale be the original scale of the data
x_posterior <- sweep(sir, MARGIN = 2, STATS = attr(X, "scaled:scale"), FUN = "*")
x_posterior <- sweep(x_posterior, MARGIN = 2, STATS = attr(X, "scaled:center"),
                     FUN = "+")
x_proposal <- sweep(x_unif, MARGIN = 2, STATS = attr(X, "scaled:scale"),
                    FUN = "*")
x_proposal <- sweep(x_proposal, MARGIN = 2,
                    STATS = attr(X, "scaled:center"), FUN = "+")

map_vars <- c("mtco" = "MTCO (°C)", "gdd5" = "GDD5 (°C days)",
              "aet.pet" = "AET/PET")

# Create data frame
list_dat <- lapply(seq_len(dim(x_posterior)[3]), function(i) x_posterior[,,i])
data_posterior <- do.call(rbind, list_dat)
dim(data_posterior)
colnames(data_posterior) <- chosen_covariates
data_posterior <- cbind(age = rep(age, each = dim(x_posterior)[1L]), data_posterior)
data_posterior <- data.frame(data_posterior)
head(data_posterior)

data_proposal <- data.frame(x_proposal)
colnames(data_proposal) <- chosen_covariates

# Plot for given years
age_t <-  9.8309
age_t_1 <- 23.9050 #27.4350
chosen_id_t <- which(age >= age_t - 0.01 & age <= age_t + 0.01)[1L]
chosen_id_t_1 <- which(age >= age_t_1 - 0.01 & age <= age_t_1 + 0.01)[1L]
chosen_age_t <- age[chosen_id_t]
chosen_age_t_1 <- age[chosen_id_t_1]

plot_bivariate_climates <- function(data_posterior, data_proposal, caption) {

  medians <- apply(data_posterior[, c("gdd5", "mtco", "aet.pet")], 2, median)

  p1 <- ggplot(data_posterior, aes(x = gdd5, y = mtco)) +
    geom_density2d(data = data_proposal, col = "brown", alpha = 0.4) +
    geom_density2d(col = "black", alpha = 0.6) +
    geom_point(size = 0.1, col = scales::alpha("black", 0.4)) +
    annotate("point", x = medians["gdd5"], y = medians["mtco"], colour = "blue",
             shape = 4, size = 4) +
    scale_x_continuous(breaks = scales::pretty_breaks(6)) +
    scale_y_continuous(breaks = scales::pretty_breaks(6)) +
    labs(x = "GDD5 (°C days)", y = "MTCO (°C)")

  p2 <- ggplot(data_posterior, aes(x = gdd5, y = aet.pet)) +
    geom_density2d(data = data_proposal, col = "brown", alpha = 0.4) +
    geom_density2d(col = "black", alpha = 0.6) +
    geom_point(size = 0.1, col = scales::alpha("black", 0.4)) +
    annotate("point", x = medians["gdd5"], y = medians["aet.pet"], colour = "blue",
             shape = 4, size = 4) +
    scale_x_continuous(breaks = scales::pretty_breaks(6)) +
    scale_y_continuous(breaks = scales::pretty_breaks(6)) +
    labs(x = "GDD5 (°C days)", y = "AET/PET")

  p3 <- ggplot(data_posterior, aes(x = mtco, y = aet.pet)) +
    geom_density2d(data = data_proposal, col = "brown", alpha = 0.4) +
    geom_density2d(col = "black", alpha = 0.6) +
    geom_point(size = 0.1, col = scales::alpha("black", 0.4)) +
    annotate("point", x = medians["mtco"], y = medians["aet.pet"], colour = "blue",
             shape = 4, size = 4) +
    scale_x_continuous(breaks = scales::pretty_breaks(6)) +
    scale_y_continuous(breaks = scales::pretty_breaks(6)) +
    labs(x = "MTCO (°C)", y = "AET/PET")

  # Grid
  plot_row <- cowplot::plot_grid(p1, p2, p3, nrow = 1, align = "hv")

  # Add the title
  title <- ggdraw() +
    draw_label(caption, fontface = 'bold', x = .45, hjust = 0) +
    theme(plot.margin = margin(0, 0, 0, 7))

  plot_grid(title, plot_row, ncol = 1, rel_heights = c(0.05, 1))
}
p_age_t_1 <- plot_bivariate_climates(data_posterior = dplyr::filter(data_posterior,
                                                                    age == chosen_age_t_1),
                                     data_proposal = data_proposal,
                                     caption = sprintf("%.3f ka BP", chosen_age_t_1))
p_age_t <- plot_bivariate_climates(data_posterior = dplyr::filter(data_posterior,
                                                                  age == chosen_age_t),
                                   data_proposal = data_proposal,
                                   caption = sprintf("%.3f ka BP", chosen_age_t))

p_final <- cowplot::plot_grid(p_age_t_1, p_age_t, ncol = 1, align = "hv",
                              labels = "AUTO")

save_plot(filename = file.path(results_dir, "kde2d_reconstructions_ages_9_23.pdf"),
          plot = p_final, base_height = 8.5)
save_plot(filename = file.path(results_dir, "kde2d_reconstructions_ages_9_23.png"),
          plot = p_final, base_height = 8.5, bg = "white")

# Plot time series -------------------------------------------------------------

# Organise data
list_data <- vector(mode = "list", length = p)
for (k in seq_len(p)) {
  post <- x_posterior_rescale[, k, ]
  list_data[[k]] <- data.frame(age = age,
                               variable = chosen_covariates[k],
                               mean = colMeans(post),
                               median = apply(post, 2, median),
                               lower_95 = apply(post, 2, quantile, probs = 0.025),
                               upper_95 = apply(post, 2, quantile, probs = 0.975),
                               lower_50 = apply(post, 2, quantile, probs = 0.75),
                               upper_50 = apply(post, 2, quantile, probs = 0.25)
                               )
}
data_posterior <- do.call(rbind, list_data)
data_posterior$variable <- map_vars[data_posterior$variable]
data_posterior$variable <- forcats::fct_relevel(data_posterior$variable,
                                                "GDD5 (°C days)", "MTCkde2d_reconstruction_age_9_23O (°C)")
head(data_posterior)
data_posterior <- dplyr::as_tibble(data_posterior)
dplyr::filter(data_posterior, age >= 15.2, age <= 18.1, variable == "AET/PET")


# Plot with ggplot
p_ts <- ggplot(data_posterior, aes(x = age, y = median)) +
  facet_wrap(~variable, scales = "free_y", ncol = 1L) +
  geom_ribbon(aes(ymin = lower_95, ymax = upper_95), fill = "grey80", alpha = 0.6) +
  # geom_ribbon(aes(ymin = lower_50, ymax = upper_50), fill = "grey80", alpha = 0.8) +
  geom_line(linewidth = .5) +
  scale_x_continuous(breaks = scales::pretty_breaks(20)) +
  scale_y_continuous(breaks = scales::pretty_breaks(6)) +
  labs(x = "Age (ka BP)", y = "")
# seq(0, 130, by = 10)
save_plot(filename = file.path(results_dir, "reconstructions.pdf"), plot = p_ts,
          base_height = 8)
save_plot(filename = file.path(results_dir, "reconstructions.png"), plot = p_ts,
          base_height = 8, bg = "white")
