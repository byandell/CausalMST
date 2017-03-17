ParametricIUCMST <- function(Z) {
  pv <- matrix(NA, 4, 4)
  for (i in 1 : 3) {
    for(j in (i + 1) : 4) {
      pv[i, j] <- pnorm(Z[i, j], lower.tail = FALSE)
    }
  }
  pval.1 <- max(pv[1, 2], pv[1, 3], pv[1, 4])
  pval.2 <- max(1 - pv[1, 2], pv[2, 3], pv[2, 4])
  pval.3 <- max(1 - pv[1, 3], 1 - pv[2, 3], pv[3, 4])
  pval.4 <- max(1 - pv[1, 4], 1 - pv[2, 4], 1 - pv[3, 4])
  c(pval.1, pval.2, pval.3, pval.4)
}
