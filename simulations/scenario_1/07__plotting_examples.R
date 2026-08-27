rm(list = ls())
library(ggplot2)
library(cowplot)
theme_set(theme_cowplot() + background_grid() +
            theme(legend.position = "top"))

d <- 4
r <- 1L

# Path
path_local <- sprintf("./simulations/scenario_1/d=%i", d)
path_data <- file.path(path_local, "data")
path_results <- file.path(path_local, "results",  sprintf("replica=%i", r))
path_zanim_ln_bart <- file.path(path_results, "zanim_ln_bart")
path_dm_bummer <- file.path(path_results, "dm_bummer")
path_dm_gp <- file.path(path_results, "dm_gp")
  
# Import data
list_data <- readRDS(file.path(path_data, sprintf("data_replica=%i.rds", r)))
id_test <- list_data$id_test
n_test <- length(id_test)
Y_test <- list_data$Y[id_test, ]
X_test <- list_data$X[id_test, , drop = FALSE]
rangeX <- range(list_data$X)
data_sim <- list_data$df[!(list_data$df$id %in% id_test), ]
data_sim$id <- rep(seq_len(nrow(list_data$Y) - n_test), each = d)

# Load inverse posterior
x_sir <- readRDS(file.path(path_zanim_ln_bart, "sir.rds"))
x_ess <- readRDS(file.path(path_zanim_ln_bart, "ess.rds"))
x_cess <- readRDS(file.path(path_zanim_ln_bart, "cess.rds"))
x_dm_bum <- readRDS(file.path(path_dm_bummer, "ip_dm_bummer.rds"))
x_dm_gp <- readRDS(file.path(path_dm_gp, "ip_dm_gp.rds"))

# Load posterior of theta
data_theta_zanim_ln <- readRDS(file.path(path_zanim_ln_bart, "posterior_theta.rds"))
data_theta_dm_bum <- readRDS(file.path(path_dm_bummer, "posterior_theta.rds"))
data_theta_dm_gp <- readRDS(file.path(path_dm_gp, "posterior_theta.rds"))
data_theta_zanim_ln$method <- "ZANIM-LN-BART"
data_theta_dm_bum$method <- "DM-BUMMER"
data_theta_dm_gp$method <- "DM-GP"
data_post_theta <- rbind(data_theta_zanim_ln, data_theta_dm_bum, data_theta_dm_gp)

# Choose specific cases
NDPOST <- 2000L
chosen_ids <- c(7L, 2L, 3L, 65L)
label_ids <- apply(Y_test[chosen_ids,], 1,
                   function(x) paste0("bold(y)[s] == (", paste0(x, collapse = ", "), ")"))
label_ids <- apply(Y_test[chosen_ids, ], 1, function(x) {
  sprintf("y[s] == '%s'",
          paste0("(", paste(x, collapse = ", "), ")"))
})
label_ids[1] <- paste0("'no zero: ' * ", label_ids[1])
label_ids[2] <- paste0("'one zero: ' * ", label_ids[2])
label_ids[3] <- paste0("'two zeros: ' * ", label_ids[3])
label_ids[4] <- paste0("'N-inflation: ' * ", label_ids[4])
# Organise data
data_x_true <- data.frame(x = X_test[chosen_ids, ], ids = chosen_ids,
                          label_ids = label_ids)
data_sir <- data.frame(x = c(x_sir[1:NDPOST,1,chosen_ids]),
                       ids = rep(chosen_ids, each = NDPOST),
                       label_ids = rep(label_ids, each = NDPOST),
                       method = "ZANIM-LN-BART (SIR)")
data_ess <- data.frame(x = c(x_ess[1:NDPOST,1,chosen_ids]),
                       ids = rep(chosen_ids, each = NDPOST),
                       label_ids = rep(label_ids, each = NDPOST),
                       method = "ZANIM-LN-BART (ESS)")
data_cess <- data.frame(x = c(x_cess[1:NDPOST,1,chosen_ids]),
                       ids = rep(chosen_ids, each = NDPOST),
                       label_ids = rep(label_ids, each = NDPOST),
                       method = "ZANIM-LN-BART (cESS)")
data_dm_gp <- data.frame(x = c(x_dm_gp[1:NDPOST,1,chosen_ids]),
                         ids = rep(chosen_ids, each = NDPOST),
                         label_ids = rep(label_ids, each = NDPOST),
                         method = "DM-GP (ESS)")
data_dm_bum <- data.frame(x = c(x_dm_bum[1:NDPOST,1,chosen_ids]),
                         ids = rep(chosen_ids, each = NDPOST),
                         label_ids = rep(label_ids, each = NDPOST),
                         method = "DM-BUMMER (ESS)")
data_post <- rbind(data_sir, data_ess, data_cess, data_dm_gp, data_dm_bum)
data_post$method <- forcats::fct_relevel(data_post$method, "ZANIM-LN-BART (SIR)", 
                                         "ZANIM-LN-BART (ESS)",
                                         "ZANIM-LN-BART (cESS)", "DM-BUMMER (ESS)")
data_post$label_ids <- forcats::fct_relevel(data_post$label_ids, label_ids)
data_x_true$label_ids <- forcats::fct_relevel(data_x_true$label_ids, label_ids)

COLORS <- c("#F8766D", "#A3A500", "#00BF7D", "#00B0F6", "#E76BF3")
# Plot


p_inverse <- ggplot(data_post) +
  facet_wrap(~label_ids, scales = "free", labeller = label_parsed) +
  stat_density(aes(x = x, col = method), geom = "line", position = "identity") +
  # geom_density(aes(x = x, col = method), show.legend = FALSE) +
  geom_point(data = data_x_true, aes(x = x, y = 0.00001)) +
  labs(x = expression(x[s]), y = "density", col = "") +
  guides(colour = guide_legend(nrow = 2, byrow = TRUE)) +
  scale_x_continuous(breaks = scales::pretty_breaks(6))
save_plot(filename = file.path(path_results, "examples_reconstruction__scenario_1.pdf"),
          plot = p_inverse, base_height = 8)

# Plot forward model fit
data_post_theta$method <- forcats::fct_relevel(data_post_theta$method, "ZANIM-LN-BART", 
                                               "DM-BUMMER")
data_post_theta$category_lab <- paste0("j == ",data_post_theta$category)
data_sim$category_lab <- paste0("j == ",data_sim$category)
p_forward <- ggplot(data = data_sim) +
  geom_line(mapping = aes(x = x, y = theta, col = "Truth", fill = "Truth"),
            linewidth = 0.8) +
  facet_wrap(~category_lab, scales = "free_y", labeller = label_parsed) +
  geom_rug(data = dplyr::filter(data_sim, total == 0L),
           mapping = aes(y = NA_real_, x = x)) +
  geom_line(data = data_post_theta, mapping = aes(x = x, y = median, col = method)) +
  geom_ribbon(data = data_post_theta,
              aes(x = x, ymin = ci_lower, ymax = ci_upper, fill = method),
              alpha = 0.3) +
  scale_color_manual(breaks = c("Truth", "ZANIM-LN-BART", "DM-BUMMER", "DM-GP"),
                     values = c("Truth" = "black",
                                "ZANIM-LN-BART" = COLORS[3L],
                                "DM-BUMMER" = COLORS[4L],
                                "DM-GP" = COLORS[5L]
                     )) +
  scale_fill_manual(breaks = c("Truth", "ZANIM-LN-BART", "DM-BUMMER", "DM-GP"),
                    values = c("Truth" = "black",
                               "ZANIM-LN-BART" = COLORS[3L],
                               "DM-BUMMER" = COLORS[4L],
                               "DM-GP" = COLORS[5L])) +
  labs(x = expression(x[s]), y = "Compositional probabilities", col = "", fill = "") +
  scale_x_continuous(breaks = scales::pretty_breaks(7))
save_plot(filename = file.path(path_results, "forward_fit__scenario_1.pdf"),
          plot = p_forward, base_height = 8)


# Grid
p_grid <- plot_grid(p_forward, p_inverse, ncol = 1, labels = "AUTO", 
                    align = 'v', axis = 'b')
p_grid
save_plot(filename = file.path(path_results, "forward_reconstruction__scenario_1_taller.pdf"),
          plot = p_grid, base_height = 10, base_width = 10)

