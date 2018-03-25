#' @export
model_chain <- function(fit1, fit2) {
  list(logLik = fit1$logLik + fit2$logLik,
       df = fit1$df + fit2$df,
       ind_logLik = fit1$ind_logLik + fit2$ind_logLik)
}