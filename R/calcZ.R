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

  LRt <- outer(ICs, ICs, function(x,y) y-x)
  LRt <- LRt[lower.tri(LRt)]
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
  names(out) <- c("ref", "alt")
  out <- dplyr::mutate(out,
                       ref = as.character(ref),
                       alt = as.character(alt))
  out$Z <- Z
  out2 <- out[,c(2,1,3)]
  out2$Z <- -Z
  names(out2) <- names(out)
  dplyr::bind_rows(out, out2)
}
