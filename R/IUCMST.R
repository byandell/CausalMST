#' @export
IUCMST <- function(models, ...) {
  result <-
    dplyr::inner_join(
      dplyr::rename(
        normIUCMST(models),
        alt.norm = alt,
        norm = pv),
      dplyr::rename(
        normJointIUCMST(models),
        alt.joint = alt,
        normJoint = pv),
      by = "ref")
  dplyr::inner_join(result,
    dplyr::rename(
        binomIUCMST(models),
        alt.binom = alt,
        binom = pv),
    by = "ref")
}
