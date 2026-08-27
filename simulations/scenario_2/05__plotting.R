rm(list = ls())
library(zanicc)
library(ggplot2)
library(cowplot)
theme_set(
  theme_cowplot(font_size = 16, font_family = "Palatino") +
    background_grid() 
)

d <- 28L
r <- 2L
p <- 2L

# Path
path_local <- sprintf("./simulations/scenario_2/d=%i_p=%i", d, p)
path_data <- file.path(path_local, "data")
path_results <- file.path(path_local, "results", sprintf("replica=%i", r))

# Import data
list_data <- readRDS(file.path(path_data, sprintf("data_replica=%i.rds", r)))
id_test <- list_data$id_test
n_test <- length(id_test)
Y_test <- list_data$Y[id_test, ]
X_test <- list_data$X[id_test, , drop = FALSE]
X_train <- list_data$X[-id_test, , drop = FALSE]

# Import posterior distribution
x_sir <- readRDS(file.path(path_results, "sir.rds"))
x_ess <- readRDS(file.path(path_results, "ess.rds"))
x_cess <- readRDS(file.path(path_results, "cess.rds"))
x_proposal <- readRDS(file.path(path_results, "x_proposal.rds"))

# Choosing some cases
zero_per_row <- apply(Y_test, 1, function(x) sum(x==0))
any(zero_per_row==0)
zero_per_row[which.min(zero_per_row)]
zero_per_row[which.max(zero_per_row)]

chosen_ids <- c(13, 70) #66 70
NDPOST <- dim(x_sir)[1L]

# Compute KDE for the convex hull of proposal = data
x1range <- range(x_proposal[, 1])
x2range <- range(x_proposal[, 2])
dens_prior <- MASS::kde2d(x_proposal[, 1], x_proposal[, 2])
# dens_prior <- MASS::kde2d(X_train[, 1], X_train[, 2])

# contour(dens_prior$x, dens_prior$y, dens_prior$z)
# contour(dens_prior2$x, dens_prior2$y, dens_prior2$z)


# Plot chosen cases
pdf(file.path(path_results, "contour_kde2d_scenario_2.pdf"), width = 7, height = 4)
par(mfrow = c(2, 3), mar = c(4, 4, 1, 1))
for (i in chosen_ids) {
  x_true <- X_test[i, ]
  # kernel densities
  dens_sir <- MASS::kde2d(x_sir[, 1,i], x_sir[,2,i], n = 100)
  dens_ess <- MASS::kde2d(x_ess[, 1,i], x_ess[,2,i], n = 100)
  dens_cess <- MASS::kde2d(x_cess[, 1,i], x_cess[,2,i], n = 100)
  
  # Get the posterior medians
  medians_sir <- apply(x_sir[, ,i], 2, median)
  medians_ess <- apply(x_ess[, ,i], 2, median)
  medians_cess <- apply(x_cess[, ,i], 2, median)
  
  # noise <- stats::rnorm(NDPOST, sd=0.01)
  
  contour(dens_prior$x, dens_prior$y, dens_prior$z,
          col = scales::alpha("brown", 0.4), main = "SIR")
  contour(dens_sir$x, dens_sir$y, dens_sir$z, add = TRUE)
  points(x_sir[, 1,i], x_sir[,2,i], col = scales::alpha("black", 0.4), cex = 0.1)
  points(x_true[1], x_true[2], col = "orange", pch = 4, cex = 2)
  points(medians_sir[1L], medians_sir[2L], col = "blue", pch = 4, cex = 2)
  mtext(expression(x[s1]), side = 1, line = 2)
  mtext(expression(x[s2]), side = 2, line = 2)
  abline(v = x_true[1], h = x_true[2], lty = 2)
  
  contour(dens_prior$x, dens_prior$y, dens_prior$z,
          col = scales::alpha("brown", 0.4), main = "ESS",
          xlim = range(x1range, x_ess[,1,i]), ylim = range(x2range, x_ess[,2,i]) )
  contour(dens_ess$x, dens_ess$y, dens_ess$z, add = TRUE)
  points(x_ess[, 1,i], x_ess[,2,i], col = scales::alpha("black", 0.4), cex = 0.1)
  points(x_true[1], x_true[2], col = "orange", pch = 4, cex = 2)
  points(medians_ess[1L], medians_ess[2L], col = "blue", pch = 4, cex = 2)
  mtext(expression(x[s1]), side = 1, line = 2)
  mtext(expression(x[s2]), side = 2, line = 2)
  abline(v = x_true[1], h = x_true[2], lty = 2)
  
  contour(dens_prior$x, dens_prior$y, dens_prior$z,
          col = scales::alpha("brown", 0.4), main = "cESS",
          xlim = range(x1range, x_cess[,1,i]), ylim = range(x2range, x_cess[,2,i]) )
  contour(dens_cess$x, dens_cess$y, dens_cess$z, add = TRUE)
  points(x_cess[, 1,i], x_cess[,2,i], col = scales::alpha("black", 0.4), cex = 0.1)
  points(x_true[1], x_true[2], col = "orange", pch = 4, cex = 2)
  points(medians_cess[1L], medians_cess[2L], col = "blue", pch = 4, cex = 2)
  mtext(expression(x[s1]), side = 1, line = 2)
  mtext(expression(x[s2]), side = 2, line = 2)
  abline(v = x_true[1], h = x_true[2], lty = 2)
  
}
graphics.off()


# Plot all
pdf(file.path(path_results, "contour_all.pdf"), width = 6, height = 3)
for (i in seq_len(n_test)) {
  cat(i, "\n")
  x_true <- X_test[i, ]
  
  # kernel densities
  dens_sir <- MASS::kde2d(x_sir[, 1,i], x_sir[,2,i], n = 100)
  dens_ess <- MASS::kde2d(x_ess[, 1,i], x_ess[,2,i], n = 100)
  dens_cess <- MASS::kde2d(x_cess[, 1,i], x_cess[,2,i], n = 100)
  
  par(mfrow = c(1, 3), mar = c(3, 3, 1, 1))
  contour(dens_prior$x, dens_prior$y, dens_prior$z,
          col = scales::alpha("brown", 0.4), main = "SIR")
  contour(dens_sir$x, dens_sir$y, dens_sir$z, add = TRUE)
  points(x_sir[, 1,i], x_sir[,2,i], col = scales::alpha("black", 0.4), cex = 0.1)
  points(x_true[1], x_true[2], col = "blue", pch = 4, cex = 2)
  abline(v = x_true[1], h = x_true[2])
  
  contour(dens_prior$x, dens_prior$y, dens_prior$z,
          col = scales::alpha("brown", 0.4), main = "ESS",
          xlim = range(x1range, x_ess[,1,i]), ylim = range(x2range, x_ess[,2,i]) )
  contour(dens_ess$x, dens_ess$y, dens_ess$z, add = TRUE)
  points(x_ess[, 1,i], x_ess[,2,i], col = scales::alpha("black", 0.4), cex = 0.1)
  points(x_true[1], x_true[2], col = "blue", pch = 4, cex = 2)
  abline(v = x_true[1], h = x_true[2])
  
  contour(dens_prior$x, dens_prior$y, dens_prior$z,
          col = scales::alpha("brown", 0.4), main = "cESS",
          xlim = range(x1range, x_cess[,1,i]), ylim = range(x2range, x_cess[,2,i]) )
  contour(dens_cess$x, dens_cess$y, dens_cess$z, add = TRUE)
  points(x_cess[, 1,i], x_cess[,2,i], col = scales::alpha("black", 0.4), cex = 0.1)
  points(x_true[1], x_true[2], col = "blue", pch = 4, cex = 2)
  abline(v = x_true[1], h = x_true[2])
  
}
graphics.off()


# ggplot style -----------------------------------------------------------------

chosen_ids <- c(13, 70) #66 70
NDPOST <- dim(x_sir)[1L]
n_sample <- dim(x_sir)[3]

# Organising data
data_sir <- do.call(rbind, lapply(seq_len(dim(x_sir)[3]), function(i) x_sir[,,i]))
data_ess <- do.call(rbind, lapply(seq_len(dim(x_ess)[3]), function(i) x_ess[,,i]))
data_cess <- do.call(rbind, lapply(seq_len(dim(x_cess)[3]), function(i) x_cess[,,i]))
colnames(data_sir) <- colnames(data_ess) <- colnames(data_cess) <- c("x1", "x2")
data_sir <- data.frame(data_sir)
data_ess <- data.frame(data_ess)
data_cess <- data.frame(data_cess)
data_sir$ids <- data_ess$ids <- data_cess$ids <- rep(seq_len(n_sample), each = NDPOST)
data_sir$method <- "SIR" 
data_ess$method <- "ESS"
data_cess$method <- "cESS"
data_posterior <- rbind(data_sir, data_ess, data_cess)
data_posterior$method <- forcats::fct_relevel(data_posterior$method, c("SIR", "ESS"))
data_proposal <- data.frame(x_proposal)
colnames(data_proposal) <- c("x1", "x2")
head(data_posterior)
head(data_proposal)



# 
data_cur_1 <- dplyr::filter(data_posterior, ids == chosen_ids[1L])
data_cur_2 <- dplyr::filter(data_posterior, ids == chosen_ids[2L])
data_true <- data.frame(X_test[chosen_ids, ])
data_medians <- data_posterior |> 
  dplyr::filter(ids %in% chosen_ids) |> 
  dplyr::group_by(ids, method) |> 
  dplyr::summarise(x1 = median(x1), x2 = median(x2), .groups = "drop")
p1 <- ggplot(data_cur_1) +
  facet_wrap(~method) +
  geom_density2d(data = data_proposal, aes(x = x1, y = x2),
                 col = "brown", alpha = 0.4) +
  geom_point(aes(x = x1, y = x2), size = 0.1, col = scales::alpha("black", 0.4)) +
  geom_density2d(aes(x = x1, y = x2), n = 100, contour_var = "ndensity",
                 col = "black", alpha = 0.6) +
  geom_vline(xintercept = data_true$X1[1L], linetype = "dashed", alpha = 0.5) +
  geom_hline(yintercept = data_true$X2[1L], linetype = "dashed", alpha = 0.5) +
  geom_point(data = dplyr::filter(data_medians, ids == chosen_ids[1L]),
             aes(x = x1, y = x2), colour = "blue", shape = 4, size = 4) +
  annotate("point", x = data_true$X1[1L], y = data_true$X2[1L], colour = "orange",
           shape = 4, size = 4) +
  scale_x_continuous(breaks = scales::pretty_breaks(6)) +
  scale_y_continuous(breaks = scales::pretty_breaks(6)) +
  labs(x = expression(x[s1]), y = expression(x[s2]))
p2 <- ggplot(data_cur_2) +
  facet_wrap(~method) +
  geom_density2d(data = data_proposal, aes(x = x1, y = x2),
                 col = "brown", alpha = 0.4) +
  geom_point(aes(x = x1, y = x2), size = 0.1, col = scales::alpha("black", 0.4)) +
  geom_density2d(aes(x = x1, y = x2), n = 100, contour_var = "ndensity",
                 col = "black", alpha = 0.6) +
  geom_vline(xintercept = data_true$X1[2L], linetype = "dashed", alpha = 0.5) +
  geom_hline(yintercept = data_true$X2[2L], linetype = "dashed", alpha = 0.5) +
  geom_point(data = dplyr::filter(data_medians, ids == chosen_ids[2L]),
             aes(x = x1, y = x2), colour = "blue", shape = 4, size = 4) +
  annotate("point", x = data_true$X1[2L], y = data_true$X2[2L], colour = "orange",
           shape = 4, size = 4) +
  scale_x_continuous(breaks = scales::pretty_breaks(6)) +
  scale_y_continuous(breaks = scales::pretty_breaks(6)) +
  labs(x = expression(x[s1]), y = expression(x[s2]))

p_g <- plot_grid(p1, p2, align = "hv", labels = "AUTO", ncol = 1)
p_g
save_plot(filename = file.path(path_results, "ggplot_contour_kde2d_scenario.pdf"),
          plot = p_g, base_height = 8.5)
