#' Compute posterior risk metrics
#' @param x matrix with observed values.
#' @param draws array with the posterior draws.
#' @param probs vector with the probabilities for the HPD region.
posterior_risk <- function(x, draws, probs = c(0.95, 0.90, 0.80, 0.50),
                           joint_coverage = TRUE) {
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

      den <- tryCatch(ks::kde(x = post), error = function(e) NULL)
      if (is.null(den)) {
        H <- ks::Hpi(x = post, binned = FALSE, pilot = "dscalar")
        den <- ks::kde(x = post, H = H, binned = TRUE)
      }

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
    # Energy-score
    es <- scoringRules::es_sample(y = x_cur, dat = t(post))

    c(mae = sum(abs(x_cur - md)), msep = sum((x_cur - mu)^2), dmode = dmode,
      es = es, crps = mean(crps), coverage, hdis_95, hdis_50, crps)
  })
  rowMeans(do.call(cbind, l))
}
