#' Log likelihood calculations
#' 
#' Log likelihood calculations used in CMST test construction.
#' @param y outcome vector
#' @param X design matrix
#' @param ... other parameters possibly used
#' 
#' @return list of
#' \itemize{
#' \item{logLik} log likelihood
#' \item{ind_logLik} vector of individual log likelihood components (density at y ~ Xb)
#' \item{d} model degrees of freedom
#' \item{RSS} residual sums of squares (optional)
#' \
#' }
logLik_calcs <- function(y, X, ...) {
  n <- length(y)
  
  dX <- ncol(X)
  qrX <- qr(X)
  b <- qr.coef(qrX, y)
  RSS <- y - X %*% b
  RSS <- crossprod(RSS, RSS)
  logLik <- as.vector(- (n/2) - (n/2) * log(2 * pi) - (n/2) * log(RSS/n))
  ss <- RSS/n
  ind_logLik <- dnorm(y, X %*% b, sqrt(ss), log = TRUE)
  list(logLik = logLik, ind_logLik = ind_logLik, df = dX, RSS = RSS)
}
