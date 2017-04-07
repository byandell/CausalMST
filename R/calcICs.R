#' @export
#'
calcICs <- function(models, flavor = c("B","A")) {

  flavor <- match.arg(flavor)
  n_ind <- nrow(models$indLR)
  penalty <- switch(flavor,
                    A = 2,
                    B = log(n_ind))

  -2 * models$LR + models$df * penalty
}
calcBICs <- function(n, model.dim, loglik) {
  calcIDs(n, model.dim, loglik, "B")
}
calcAICs <- function(n, model.dim, loglik) {
  calcIDs(n, model.dim, loglik, "A")
}
