#' @export
#' 
ParametricIUCMST <- function(models, Zscores = calcZ(models)) {
  
  # Compare each model with all others and get max pvalue.
  # Need to handle left and right in opposite ways.
  ref <- unique(Zscores$ref)
  Zscores <-
    dplyr::ungroup(
      dplyr::summarize(
        dplyr::group_by(Zscores, ref),
        alt = alt[which.max(pv)][1],
        pv = max(pv)))

  out <- Zscores$pv
  names(out) <- Zscores$ref
  out[ref]
}
