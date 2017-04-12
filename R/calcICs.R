#' @export
#'
calcICs <- function(models, flavor = c("B","A")) {

  flavor <- match.arg(flavor)
  n_ind <- nrow(models$indLR)
  -2 * models$LR + models$df * penalty(n_ind, flavor)
}
calcBICs <- function(n, model.dim, loglik) {
  calcICs(n, model.dim, loglik, "B")
}
calcAICs <- function(n, model.dim, loglik) {
  calcICs(n, model.dim, loglik, "A")
}
penalty <- function(n_ind, flavor) {
  switch(flavor,
         A = 2,
         B = log(n_ind))
}
