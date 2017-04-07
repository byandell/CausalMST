#' @export
model_chain <- function(fit1, fit2) {
  list(loglik = fit1$log.lik + fit2$log.lik,
       df = fit1$df + fit2$df,
       vec.logLik = fit1$vec.log.lik + fit2$vec.log.lik)
}