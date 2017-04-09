#' @export
#' 
ParametricIUCMST <- function(models, Zscores = calcZ(models)) {
  
  Zscores$pv <- pnorm(Zscores$Z, lower.tail = FALSE)

  # Compare each model with all others and get max pvalue.
  # Need to handle left and right in opposite ways.
  Zscores <-
    dplyr::ungroup(
      dplyr::summarize(
        dplyr::group_by(Zscores, ref),
        alt = alt[which.max(pv)][1],
        pv = max(pv)))

  out <- Zscores$pv
  names(out) <- Zscores$ref
  out
}
