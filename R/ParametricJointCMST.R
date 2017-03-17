ParametricJointCMST <- function(Z, Cor.hat) {
  z <- min(Z[1, 2], Z[1, 3], Z[1, 4])
  pval.1 <- 1 - pmnorm(c(z, z, z), c(0, 0, 0), Cor.hat[[1]])
  z <- min(- Z[1, 2], Z[2, 3], Z[2, 4])
  pval.2 <- 1 - pmnorm(c(z, z, z), c(0, 0, 0), Cor.hat[[2]])
  z <- min(-Z[1, 3], -Z[2, 3], Z[3, 4])
  pval.3 <- 1 - pmnorm(c(z, z, z), c(0, 0, 0), Cor.hat[[3]])
  z <- min(-Z[1, 4], -Z[2, 4], -Z[3, 4])
  pval.4 <- 1 - pmnorm(c(z, z, z), c(0, 0, 0), Cor.hat[[4]])
  c(pval.1, pval.2, pval.3, pval.4)
}
