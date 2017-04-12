corHat <- function(S.hat) {
  tmp <- corpcor::is.positive.definite(S.hat)
  if (!tmp) {
    S.hat <- corpcor::make.positive.definite(S.hat)
  }
  cov2cor(S.hat)
}
