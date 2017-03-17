calcBICs <- function(n, model.dim, loglik) {
  -2 * loglik + model.dim * log(n)
}
calcAICs <- function(n, model.dim, loglik) {
  -2 * loglik + 2 * model.dim
}
