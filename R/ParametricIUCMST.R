#' @export
#' 
ParametricIUCMST <- function(object) {
  
  object$pv <- pnorm(object$Z, lower.tail = FALSE)

  left <- dplyr::rename(
    dplyr::ungroup(
      dplyr::summarize(
        dplyr::group_by(object, left),
        right = right[which.max(pv)][1],
        pvr = max(pv))),
    model = left)
  right <- dplyr::rename(
    dplyr::ungroup(
      dplyr::summarize(
        dplyr::group_by(object, right),
        left = left[which.min(pv)][1],
        pvl = max(1 - pv))),
    model = right)

  object <- dplyr::full_join(left,right, by="model")
  object <- dplyr::mutate(object,
                          pvr = ifelse(is.na(pvr), 0, pvr),
                          pvl = ifelse(is.na(pvl), 0, pvl),
                          pvalue = pmax(pvr, pvl),
                          other = ifelse(pvr > pvl, right, left))
  out <- object$pvalue
  names(out) <- object$model
  out
}
