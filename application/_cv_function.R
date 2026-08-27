#' Compute posterior risk metrics
#' @param x matrix with observed values.
#' @param draws array with the posterior draws.
#' @param probs vector with the probabilities for the HPD region.
posterior_risk <- function(x, draws, probs = c(0.95, 0.50), binned_kde = TRUE,
                           pilot_kde = "dscalar", joint_coverage = TRUE) {
  n <- nrow(x)
  p <- ncol(x)
  stopifnot(n == dim(draws)[3L])
  l <- lapply(seq_len(n), function(i) {
    cat(i, "of", n, "\n")
    x_cur <- x[i, ]
    post <- as.matrix(draws[, , i])
    mu <- colMeans(post)
    md <- apply(post, 2, median)
    # Compute marginal HDIs
    hdis_95 <- coda::HPDinterval(coda::as.mcmc(post), prob = 0.95)
    hdis_95 <- (x_cur >= hdis_95[, "lower"]) & (x_cur <= hdis_95[, "upper"])
    hdis_50 <- coda::HPDinterval(coda::as.mcmc(post), prob = 0.50)
    hdis_50 <- (x_cur >= hdis_50[, "lower"]) & (x_cur <= hdis_50[, "upper"])
    names(hdis_95) <- paste0("coverage_95_x", 1:p)
    names(hdis_50) <- paste0("coverage_50_x", 1:p)
    
    # KDE estimate
    dmode <- coverage <- NULL
    if (joint_coverage) {
      H <- tryCatch(ks::Hpi(x = post, binned = binned_kde, pilot = pilot_kde), error = function(e) NULL)
      test <- tryCatch(expr = chol(H), error = function(e) NULL)
      if (is.null(test)) H <- ks::Hpi(x = post, binned = FALSE, pilot = "dscalar")
      den <- ks::kde(x = post, H = H, binned = TRUE)
      # interpolation over posterior draws and computation of the threshold
      fdraws <- predict(den, x = post)
      k <- quantile(fdraws, prob = 1.0 - probs)
      # interpolation to obtain the density at observed values of x
      fx <- predict(den, x = x_cur)
      coverage = as.integer(fx >= k)
      names(coverage) <- paste0("coverage_", 100*probs)
      # mode
      index <- which.max(fdraws)
      mo <- post[index,]
      dmode <- sum((x_cur - mo)^2)
    }
    # mean (over covariates) crps
    crps <- sapply(seq_len(p), function(j) {
      scoringRules::crps_sample(y = x_cur[j], dat = t(post[,j]))
    })
    names(crps) <- paste0("crps_", 1:p)
    # energy score
    es <- scoringRules::es_sample(y = x_cur, dat = t(post))
    
    c(mae = sum(abs(x_cur - md)), msep = sum((x_cur - mu)^2),
      dmode = dmode, es = es, crps = mean(crps), coverage, hdis_95, hdis_50, crps)
  })
  rowMeans(do.call(cbind, l))
}


#' Perform the k-fold cross-validation for inverse posterior using ZANIM-LN-BART
#' and the SIR algorithm
#' @param Y matrix with counts.
#' @param X matrix with covariates.
#' @param folds vector indicating which fold each observation belongs.
#' @param k integer indicating which fold to run.
#' @param models_parameters list with `ndpost`, `ntrees_theta` and `ntrees_zeta`.
#' @param x_proposal matrix with the uniform proposal for generated previously.
#' @param results_dir string with folder path to save the results.
#' @param keep_draws logical. If keep the draws or not. Default is false.
cv_zanim_ln_bart <- function(Y, X, folds, k, model_parameters, x_proposal,
                             results_dir, keep_draws = FALSE) {

  # Create a folder to keep the draws for the current CV-fold
  results_dir_foldk <- file.path(results_dir, paste0("fold_", k))
  forests_dir <- file.path(results_dir_foldk, "forests")
  if (!dir.exists(forests_dir)) dir.create(forests_dir, recursive = TRUE)

  # Split data
  idx_test <- which(folds == k, arr.ind = TRUE)
  Y_train <- Y[-idx_test, , drop = FALSE]
  X_train <- X[-idx_test, , drop = FALSE]
  Y_test <- Y[idx_test, , drop = FALSE]
  X_test <- X[idx_test, , drop = FALSE]

  # Hard coding just to have some bullshit results
  # Y_test <- Y_test[1:200, ]
  # X_test <- X_test[1:200, ]
  
  # Run forward model
  if (!file.exists(file.path(results_dir_foldk, "mod.rds"))) {
    zanim_ln_bart <- zanicc(Y = Y_train, X_count = X_train, X_zi = X_train,
                            model = "zanim_ln_bart",
                            ntrees_theta = model_parameters$ntrees_theta,
                            ntrees_zeta = model_parameters$ntrees_zeta,
                            ndpost = model_parameters$ndpost,
                            nskip = model_parameters$nskip,
                            keep_draws = keep_draws,
                            save_trees = TRUE, forests_dir = forests_dir)
    save_model(object = zanim_ln_bart, model_dir = results_dir_foldk)
  }
  zanim_ln_bart <- load_model(model_dir = results_dir_foldk)

  # Run SIR
  sir <- inverse_posterior_zanimlnbart(object = zanim_ln_bart,
                                       Y = Y_test,
                                       x_proposal = x_proposal,
                                       dir_posterior_fx = forests_dir,
                                       method = "sir",
                                       ndpost = 2000L)
  saveRDS(object = sir, file = file.path(results_dir_foldk, "sir.rds"))

  # Compute and write metrics
  sir_metrics <- posterior_risk(x = X_test, draws = sir)
  write.table(x = data.frame(k, sir_metrics[1], sir_metrics[2], sir_metrics[3],
                             sir_metrics[4], sir_metrics[5], crps_sir),
              file = file.path(results_dir, "metrics_sir.txt"), append = TRUE,
              row.names = FALSE, col.names = FALSE)
  return(1L)
}
