rm(list = ls())
library(zanicc)

# Settings
NTREES_THETA <- 100L
NTREES_ZETA <- 100L
NDPOST <- 4000L
NSKIP <- 4000L

# Paths
path_local <- "./application"
results_dir <- file.path(path_local, "results", "forward_module")
forests_dir <- file.path(results_dir, "forests")
if (!dir.exists(forests_dir)) dir.create(forests_dir, recursive = TRUE)

# Load modern data
data("pollen_climate")
Y <- pollen_climate$Y
X <- pollen_climate$X[, c("gdd5", "mtco", "aet.pet")]
X <- scale(X)

cat("Run model\n\n")
# Forward module: Fit ZANIM-LN-BART model:
if (!file.exists(file.path(results_dir, "mod.rds"))) {
  zanim_ln_bart <- zanicc(Y = Y, X_count = X, X_zi = X,
                          model = "zanim_ln_bart", ntrees_theta = NTREES_THETA,
                          ntrees_zeta = NTREES_ZETA, ndpost = NDPOST,
                          nskip = NSKIP, keep_draws = FALSE,
                          save_trees = TRUE, forests_dir = forests_dir)
  save_model(object = zanim_ln_bart, model_dir = results_dir)
}
