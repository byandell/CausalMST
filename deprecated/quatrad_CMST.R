#' Causal Model Selection Tests for mediation quatrads
#'
#' @param models Object from \code{\link{mediationModels}}
#' @param test Type of CMST test to perform
#' @param threshold Threshold for secondary test of best model across quatrads
#'
#' @export
#' @importFrom dplyr bind_rows filter inner_join rename
#' @importFrom purrr map
#'
quatrad_CMST <- function(models, test = c("wilc","binom","joint","norm"),
                         threshold = 0.1) {
  test <- match.arg(test)
  testfn <- switch(test,
                   wilc = CausalMST::wilcIUCMST,
                   binom = CausalMST::binomIUCMST,
                   joint = CausalMST::normJointIUCMST,
                   norm = CausalMST::normIUCMST)

  node_id <- quatrads()
  nodes <- node_id$nodes
  node_id <- node_id$node_id

  tmpfn <- function(x, models) {
    models <- subset(models, x)
    dplyr::filter(
      testfn(models),
      pv == min(pv))
  }

  models_pv <-
    dplyr::bind_rows(
      purrr::map(node_id,
                 tmpfn, models),
      .id = "role")
  
  ref <- (dplyr::filter(models_pv, pv <= threshold))$ref
  if(length(ref) < 2) {
    dplyr::mutate(models_pv,
                  best.pv = 1,
                  best.alt = alt)
  } else {
    best <-
      dplyr::rename(
        testfn(
          subset(models,
            unique((dplyr::filter(nodes,
                                  model %in% ref))$id))),
        best.pv = pv,
        best.alt = alt)
    dplyr::left_join(models_pv, best, by = "ref")
  }
}
