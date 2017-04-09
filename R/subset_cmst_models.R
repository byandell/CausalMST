#' @export
subset.cmst_models <- function(object, model_names, ...) {
  object$LR <- object$LR[model_names]
  object$indLR <- object$indLR[model_names]
  object$df <- object$df[model_names]
  object
}