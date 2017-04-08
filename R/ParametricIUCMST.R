#' @export
#' 
ParametricIUCMST <- function(models, Zscores = calcZ(models)) {
  
  Zscores$pv <- pnorm(Zscores$Z, lower.tail = FALSE)

  # Compare each model with all others and get max pvalue.
  # Need to handle left and right in opposite ways.
  left <-
    dplyr::ungroup(
      dplyr::summarize(
        dplyr::group_by(
          dplyr::rename(Zscores, 
                        model = left),
          model),
        right = right[which.max(pv)][1],
        pvr = max(pv))) # Use upper tail.
    
  right <-
    dplyr::ungroup(
      dplyr::summarize(
        dplyr::group_by(
          dplyr::rename(Zscores,
                        model = right),
          model),
        left = left[which.min(pv)][1],
        pvl = max(1 - pv))) # Use lower tail.

  # Join left and right model comparisons; find max pvalue for each model.
  object <- 
    dplyr::mutate(
      dplyr::full_join(left, right, by="model"),
      pvr = ifelse(is.na(pvr), 0, pvr),
      pvl = ifelse(is.na(pvl), 0, pvl),
      pvalue = pmax(pvr, pvl),
      other = ifelse(pvr > pvl, right, left))
  
  out <- object$pvalue
  names(out) <- object$model
  out
}
