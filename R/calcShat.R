calcShat <- function(vec.logLik) {
  n.ind <- nrow(vec.logLik)
  vec.LR <- matrix(NA, n.ind, 6, 
                   dimnames = list(NULL, c("12", "13", "14", "23", "24", "34")))
  ip <- matrix(c(1, 1, 1, 2, 2, 3, 2, 3, 4, 3, 4, 4), 6, 2)
  for (i in 1 : 6) {
    vec.LR[, i] <- vec.logLik[, ip[i, 1]] - vec.logLik[, ip[i, 2]]
  }
  S.hat <- (1 - 1 / n.ind) * cov(vec.LR)
  dimnames(S.hat) <- list(dimnames(vec.LR)[[2]], dimnames(vec.LR)[[2]])
  S.hat
}
