#' @export
quatrad_CMST <- function(models, test = c("wilc","binom","joint","norm")) {
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
                 tmpfn, models_par$models),
      .id = "quatrad")

  ref <- (dplyr::filter(models_pv, pv <= 0.5))$ref

  if(length(ref) > 1) {
    best <-
      dplyr::rename(
        CausalMST::wilcIUCMST(
          subset(models_par$models,
            unique((dplyr::filter(nodes,
                                  model %in% ref))$id))),
        best.pv = pv,
        best.alt = alt)
  }
  dplyr::inner_join(models_pv, best, by = "ref")
}
