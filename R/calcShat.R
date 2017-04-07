#' @importFrom dplyr bind_cols
#' @importFrom purrr map
#' @export
#'
calcShat <- function(indLR) {
  
  # Might have object with indLR as element.
  if(!(is.data.frame(indLR) | is.matrix(indLR)))
    indLR <- indLR$indLR

  d <- dim(indLR)
  indLR <- as.data.frame(indLR)

  LR <-
    dplyr::bind_cols(
      purrr::map(indLR,
                 function(x,y) x-y,
                 indLR))
  names(LR) <- paste(rep(names(indLR), each = ncol(indLR)),
                     names(LR),
                     sep = ":")
  LR <- LR[, rep(seq(d[2]), each = d[2]) < rep(seq(d[2]), d[2])]

  #Shat
  (1 - 1 / d[1]) * cov(LR)
}
