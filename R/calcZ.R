calcZ <- function(S.hat, ICs, n.ind) {
  Z.ic <- matrix(NA, 4, 4)
  ii <- 1
  for (i in 1 : 3) {
    for (j in (i + 1) : 4) {
      LRt <- -0.5 * (ICs[i] - ICs[j])
      Z.ic[i, j] <- LRt / sqrt(S.hat[ii, ii] * n.ind)
      ii <- ii + 1
    }
  }
  Z.ic
}
