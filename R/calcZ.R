#' @export
#'
#' @importFrom purrr simplify_all transpose
#' @importFrom stringr str_split
#' @importFrom dplyr mutate
#'
calcZ <- function(models,
                  S.hat = calcShat(models),
                  ICs = calcICs(models, flavor),
                  flavor = c("B","A"),
                  ...) {

  n_ind <- nrow(models$indLR)

  LRt <- outer(ICs, ICs, function(x,y) y-x)
  LRt <- LRt[lower.tri(LRt)]
  LRt <- -LRt / (2 * sqrt(n_ind))
  LRt / sqrt(diag(S.hat))
}
