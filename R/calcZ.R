#' @export
#' 
#' @importFrom purrr simplify_all transpose
#' @importFrom stringr str_split
#' @importFrom dplyr mutate
#'
calcZ <- function(models,
                  S.hat = calcShat(models),
                  ICs = calcICs(models, flavor),
                  flavor = c("B","A")) {

  n_ind <- nrow(models$indLR)
  n_mod <- length(ICs)

  LRt <- outer(ICs, ICs, function(x,y) x-y)
  LRt <- LRt[rep(seq(n_mod), each = n_mod) < rep(seq(n_mod), n_mod)]
  LRt <- -LRt / (2 * sqrt(n_ind))
  Z <- LRt / sqrt(diag(S.hat))
  
  left_right(Z)
}

left_right <- function(Z) {
  out <-
    as.data.frame(
      purrr::simplify_all(
        purrr::transpose(
          stringr::str_split(names(Z), ":"))))
  names(out) <- c("left", "right")
  out <- dplyr::mutate(out,
                       left = as.character(left),
                       right = as.character(right))
  out$Z <- Z
  out
}
