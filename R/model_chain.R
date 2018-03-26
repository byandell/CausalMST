#' @export
model_chain <- function(fit1, fit2) {
  list(LR = fit1$LR + fit2$LR,
       df = fit1$df + fit2$df,
       indLR = fit1$indLR + fit2$indLR)
}