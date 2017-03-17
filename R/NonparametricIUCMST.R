NonparametricIUCMST <- function(penalty, n, k, vec.logLik) {
  # Computes the non-parametric CMST
  # k: vector of length 4 with the model dimensions
  # vec.logLik: matrix of loglik scores (n by 4)
  kM <- matrix(rep(k, each = n), n, 4)
  if (penalty == "bic") {
    vec.penal <- vec.logLik - 0.5 * kM * log(n)/n
  }
  if (penalty == "aic") {
    vec.penal <- vec.logLik - kM/n
  }
  pv <- matrix(NA, 4, 4)
  for (i in 1 : 3) {
    for (j in (i + 1) : 4) {
      vec.ratio <- vec.penal[, i] - vec.penal[, j]
      counts <- sum(vec.ratio > 0)
      nn <- n - sum(vec.ratio == 0)
      pv[i, j] <- pbinom(counts - 1, nn, 0.5, lower.tail = FALSE)
      pv[j, i] <- pbinom(counts, nn, 0.5, lower.tail = TRUE)
    }
  }
  pval.1 <- max(pv[1, 2], pv[1, 3], pv[1, 4])
  pval.2 <- max(pv[2, 1], pv[2, 3], pv[2, 4])
  pval.3 <- max(pv[3, 1], pv[3, 2], pv[3, 4])
  pval.4 <- max(pv[4, 1], pv[4, 2], pv[4, 3])
  c(pval.1, pval.2, pval.3, pval.4)
}
