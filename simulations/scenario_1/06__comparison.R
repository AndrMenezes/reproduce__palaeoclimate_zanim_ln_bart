rm(list = ls())
library(xtable)
library(zanicc)
library(ggplot2)

d <- 4
nreplica <- 6L

# Path
path_local <- sprintf("./simulations/scenario_1/d=%i", d)
path_data <- file.path(path_local, "data")

list_metrics <- vector(mode = "list", length = nreplica)
for (r in seq_len(nreplica)) {
  
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
  
  # Load estimated posterior
  x_sir <- readRDS(file.path(path_zanim_ln_bart, "sir.rds"))
  x_ess <- readRDS(file.path(path_zanim_ln_bart, "ess.rds"))
  x_cess <- readRDS(file.path(path_zanim_ln_bart, "cess.rds"))
  x_dm_bum <- readRDS(file.path(path_dm_bummer, "ip_dm_bummer.rds"))
  x_dm_gp <- readRDS(file.path(path_dm_gp, "ip_dm_gp.rds"))
  
  # Comparing against the observed value of X
  list_ip <- list(sir = x_sir, ess = x_ess, cess = x_cess,
                  dm_gp = x_dm_gp, dm_bummer = x_dm_bum)
  metrics <- lapply(list_ip, function(dr) {
    zanicc:::compute_prediction_metrics(x = X_test, draws = dr)
  })
  list_metrics[[r]]  <- do.call(rbind, metrics)
}
perf <- do.call(rbind, list_metrics)
tab <- apply(perf, 2, function(x)  {
  tapply(x, rownames(perf), mean)
})
tab[order(tab[, 4]),]
order_rows <- c("sir", "ess", "cess", "dm_gp", "dm_bummer")
order_cols <- c("crps", "msep", "mae", "coverage_95")
tab <- tab[order_rows,order_cols]
rownames(tab) <- c("ZANIM-LN-BART (SIR)", "ZANIM-LN-BART (ESS)",
                   "ZANIM-LN-BART (cESS)", "DM-GP (ESS)", "DM-BUMMER (ESS)")
print(xtable::xtable(tab, digits = 3),
      include.rownames = TRUE)



