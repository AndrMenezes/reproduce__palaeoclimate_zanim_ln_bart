rm(list = ls())
library(zanicc)
library(ggplot2)
source("./simulations/utils.R")

d <- 28L
nrep <- 6L
p <- 3L

# Path
path_local <- sprintf("./simulations/scenario_2/d=%i_p=%i", d, p)
path_data <- file.path(path_local, "data")

list_metrics <- vector(mode = "list", length = nrep)

for (r in seq_len(nrep)) {

  cat(r, "of", nrep)

  path_results <- file.path(path_local, "results", sprintf("replica=%i", r))
  X_test <- readRDS(file.path(path_data, sprintf(path_data, "data_replica=%i.rds", r)))

  # Load estimated posterior
  x_sir <- readRDS(file.path(path_results, "sir.rds"))
  x_ess <- readRDS(file.path(path_results, "ess.rds"))
  x_cess <- readRDS(file.path(path_results, "cess.rds"))

  # Compute metrics
  list_ip <- list(sir = x_sir, ess = x_ess, cess = x_cess)
  metrics <- parallel::mclapply(list_ip, function(dr) {
    posterior_risk(x = X_test, draws = dr, joint_coverage = TRUE)
  }, mc.cores = 3)

  list_metrics[[r]] <- do.call(rbind, metrics)
}
# perf <- do.call(rbind, list_metrics)
# saveRDS(object = perf, file = file.path(path_local, "metrics.rds"))


perf <- readRDS(file = file.path(path_local, "metrics.rds"))
tab <- apply(perf, 2, function(x)  {
  tapply(x, rownames(perf), mean)
})
tab

chosen_cols <- c("es", "msep", "mae",  "coverage_95")
tab_metrics <- data.frame(tab[c("sir", "ess", "cess"), chosen_cols])

# Get the time taken
list_times <- lapply(seq_len(nrep), function(r) {
  path_results <- file.path(path_local, "results", sprintf("replica=%i", r))
  x_sir <- readRDS(file.path(path_results, "sir.rds"))
  x_ess <- readRDS(file.path(path_results, "ess.rds"))
  x_cess <- readRDS(file.path(path_results, "cess.rds"))
  c(sir_sampling = attr(x_sir, "elapsed_time")[3],
    sir_predict = attr(x_sir, "elapsed_time_predict")[3],
    ess = attr(x_ess, "elapsed_time")[3],
    cess = attr(x_cess, "elapsed_time")[3])
})
elapsed_time <- colMeans(do.call(rbind, list_times))
times <- c(sir = paste0(formatC(elapsed_time[1L], digits = 1, format = "f"),
                        " (",
                        formatC(elapsed_time[2L], digits = 1, format = "f"),
                        ")"),
           ess = formatC(elapsed_time[3L], digits = 1, format = "f") 
           ,
           cess = formatC(elapsed_time[4L], digits = 1, format = "f"))

format_sec_hour <- function(x) {
  sprintf("%dh %1d min", x %/% 3600, (x %% 3600) %/% 60)
}
sec1 <- elapsed_time[1L] + elapsed_time[2L]
sec2 <- elapsed_time[3L]
sec3 <- elapsed_time[4L]
format_sec_hour(sec1)
format_sec_hour(sec2)
format_sec_hour(sec3)
sec1
sec2
sec3

times <- c(sir = format_sec_hour(elapsed_time[1L] + elapsed_time[2L]),
           ess = format_sec_hour(elapsed_time[3L]), 
           cess = format_sec_hour(elapsed_time[4L]))
tab <- cbind(tab_metrics, time = times)
rownames(tab) <- c("SIR", "ESS", "cESS")
print(xtable::xtable(tab, digits = 3), include.rownames = TRUE)


