#' Log likelihood calculations
#' 
#' Log likelihood calculations used in CMST test construction.
#' @param y outcome vector
#' @param X design matrix
#' @param ... other parameters possibly used
#' 
#' @return list of
#' \itemize{
#' \item{log.lik} log likelihood
#' \item{vec.log.lik} vector of individual log likelihood components (density at y ~ Xb)
#' \item{d} model degrees of freedom
#' \item{RSS} residual sums of squares (optional)
#' \
#' }
logLik_calcs <- function(y, X, ...) {
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
