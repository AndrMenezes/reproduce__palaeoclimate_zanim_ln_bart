rm(list = ls())
library(zanicc)

# Paths
path_local <- "./application"
path_data <- file.path(path_local, "data")
path_results <- file.path(path_local, "results", "forward_module")


# Import data ------------------------------------------------------------------
data("pollen_climate", package = "zanicc")
X <- pollen_climate$X[,  c("gdd5", "mtco", "aet.pet")]

# Create a 3d mesh/grid inside the observed points -----------------------------

min_gdd5 <- min(X[, "gdd5"])
max_gdd5 <- max(X[, "gdd5"])
min_mtco <- min(X[, "mtco"])
max_mtco <- max(X[, "mtco"])
min_aetpet <- min(X[, "aet.pet"])
max_aetpet <- max(X[, "aet.pet"])

ran_gdd5 <- max_gdd5 - min_gdd5
ran_mtco <- max_mtco - min_mtco
ran_aetpet <- max_aetpet - min_aetpet

# Transform variables
step <- 1L
scale_factor_gdd5 <- scale_factor_mtco <- 30
scale_factor_aetpet <- 20
gdd5_scale <- scale_factor_gdd5 * (X[, "gdd5"] - min_gdd5) / ran_gdd5
mtco_scale <- scale_factor_mtco * (X[, "mtco"] - min_mtco) / ran_mtco
aetpet_scale <- scale_factor_aetpet * (X[, "aet.pet"] - min_aetpet) / ran_aetpet

# Create an uniform grid on the transformed variables
X_grid <- expand.grid(gdd5 = seq(0, scale_factor_gdd5, by = step),
                      mtco = seq(0, scale_factor_mtco, by = step),
                      aet.pet = seq(0, scale_factor_aetpet, by = step))
nrow(X_grid)

keep <- logical(length = nrow(X_grid))
for (i in seq_len(nrow(X_grid))) {
  cond_gdd5 <- (gdd5_scale > X_grid[i, 1L] - step) & (gdd5_scale < X_grid[i, 1L] + step)
  cond_mtco <- (mtco_scale > X_grid[i, 2L] - step) & (mtco_scale < X_grid[i, 2L] + step)
  cond_aetpet <- (aetpet_scale > X_grid[i, 3L] - step) & (aetpet_scale < X_grid[i, 3L] + step)
  keep[i] <- any(cond_gdd5 & cond_mtco & cond_aetpet)
}
sum(keep)
mean(keep)
X_grid[, "gdd5"] <- X_grid[, "gdd5"] * ran_gdd5 / scale_factor_gdd5 + min_gdd5
X_grid[, "mtco"] <- X_grid[, "mtco"] * ran_mtco / scale_factor_mtco + min_mtco
X_grid[, "aet.pet"] <- X_grid[, "aet.pet"] * ran_aetpet / scale_factor_aetpet + min_aetpet

# Bivariate plots --------------------------------------------------------------

png(filename = file.path(path_results, "discrete_3d_grid.png"), width = 800, height = 700)
par(mfrow = c(1, 3))
plot(X[, "gdd5"], X[, "mtco"], main = "GDD5 vs MTCO")
points(X_grid[, "gdd5"], X_grid[, "mtco"], col = "blue", pch = 19)
points(X_grid[keep, "gdd5"], X_grid[keep, "mtco"], col = "red", pch = 19)

plot(X[, "gdd5"], X[, "aet.pet"], main = "GDD5 vs AET.PET")
points(X_grid[, "gdd5"], X_grid[, "aet.pet"], col = "blue", pch = 19)
points(X_grid[keep, "gdd5"], X_grid[keep, "aet.pet"], col = "red", pch = 19)

plot(X[, "aet.pet"], X[, "mtco"], main = "MTCO vs AET.PET")
points(X_grid[, "aet.pet"], X_grid[, "mtco"], col = "blue", pch = 19)
points(X_grid[keep, "aet.pet"], X_grid[keep, "mtco"], col = "red", pch = 19)
graphics.off()

head(X_grid)

# Save
saveRDS(object = X_grid[keep, ], file = file.path(path_data, "X_3d_buffer.rds"))
