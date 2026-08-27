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
Y_train <- list_data$Y[-id_test, ]
X_train <- list_data$X[-id_test, , drop = FALSE]
true_thetas <- list_data$true_thetas[-id_test, ]
true_varthetas <- list_data$true_varthetas[-id_test, ]
true_zetas <- list_data$true_zetas[-id_test, ]
d <- ncol(Y_train)

# Fit forward model
NDPOST <- 5000L
NSKIP <- 5000L
NTREES_THETA <- 100L
NTREES_ZETA <- 100L

if (!file.exists(file.path(path_results, "mod.rds"))) {

zanim_ln_bart <- zanicc(Y = Y_train, X_count = X_train, X_zi = X_train,
                        model = "zanim_ln_bart", ntrees_theta = NTREES_THETA,
                        ntrees_zeta = NTREES_ZETA, ndpost = NDPOST,
                        nskip = NSKIP, save_trees = TRUE, keep_draws = FALSE,
                        forests_dir = forests_dir)
save_model(object = zanim_ln_bart, model_dir = path_results)

}
