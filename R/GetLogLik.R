GetLogLik <- function(geno, pheno, resp, addcov, intcov) {
  X <- CreateDesignMatrix(geno, pheno, addcov, intcov)
  logLik_calcs(pheno[[resp]], X)
}

logLik_calcs <- function(y, X) {
  n <- length(y)
  
  dX <- ncol(X)
  qrX <- qr(X)
  b <- qr.coef(qrX, y)
  RSS <- crossprod(y - X %*% b, y - X %*% b)
  log.lik <- as.vector(- (n/2) - (n/2) * log(2 * pi) - (n/2) * log(RSS/n))
  ss <- RSS/n
  vec.log.lik <- dnorm(y, X %*% b, sqrt(ss), log = TRUE)
  list(log.lik = log.lik, vec.log.lik = vec.log.lik, d = dX, RSS = RSS)
}
