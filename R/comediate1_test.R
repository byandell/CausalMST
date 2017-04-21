# Comediation tests
#
#' Develop mediation models from driver, target and mediator
#'
#' @param driver vector or matrix with driver values
#' @param target vector or 1-column matrix with target values
#' @param mediator matrix with mediator values by column
#' @param fitFunction function to fit models with driver, target and mediator
#' @param kinship optional kinship matrix among individuals
#' @param cov_tar optional covariates for target
#' @param cov_med optional covariates for mediator
#' @param annotation Table with annotation, with \code{id}
#' agreeing with column names of \code{mediator}.
#' @param test Type of CMST test.
#' @param pos Position of driver.
#' @param ... additional parameters
#'
#' @importFrom stringr str_replace
#' @importFrom qtl2scan fit1 get_common_ids
#' @importFrom dplyr arrange bind_rows filter left_join rename
#'
#' @export
#'
comediate1_test <- function(driver, target, mediator, fitFunction,
                          kinship=NULL, cov_tar=NULL, cov_med=NULL,
                          annotation, test = c("wilc","binom","joint","norm"),
                          pos = NULL,
                          ...) {

  test <- match.arg(test)
  testfn <- switch(test,
                   wilc = CausalMST::wilcIUCMST,
                   binom = CausalMST::binomIUCMST,
                   joint = CausalMST::normJointIUCMST,
                   norm = CausalMST::normIUCMST)

  pos_t <- pos

  scan_max <- fitFunction(driver, target, kinship, cov_tar)
  lod_t <- scan_max$lod

  commons <- common_data(driver, target, mediator,
                         kinship, cov_tar, cov_med)
  driver <- commons$driver
  target <- commons$target
  mediator <- commons$mediator
  kinship <- commons$kinship
  cov_tar <- commons$cov_tar
  cov_med <- commons$cov_med

  cmst_fit <- function(x, driver) {
    # Force x (= mediator column) to be matrix.
    x <- as.matrix(x)
    rownames(x) <- rownames(driver)
    colnames(x) <- "M"
    # Fit mediation models.
    models_par <- mediationModels(driver, target, x,
                                  qtl2scan::fit1,
                                  kinship, cov_tar, cov_medi,
                                  common = TRUE)
    # CMST on quatrads
    out <- dplyr::filter(
      testfn(subset(models_par$models, 1:4)),
      pv == min(pv))
    # Mediation LOD
    med_lod <- sum(models_par$comps$LR[c("t.d_t", "mediation")]) / log(10)
    # Mediator LOD
    medor_lod <- models_par$comp$LR["m.d_m"] / log(10)
    out$mediation <- med_lod
    out$mediator <- medor_lod

    out
  }

  best <- vector("list", ncol(mediator))
  names(best) <- colnames(mediator)
  for(i in names(best)) {
    covi <- unlist(dplyr::filter(annotation, pheno == i)[, colnames(cov_med)])
    cov_medi <- cov_med[, covi, drop = FALSE]
    best[[i]] <- cmst_fit(mediator[,i], driver)
  }

  best <- dplyr::rename(
    dplyr::bind_rows(best, .id = "pheno"),
    triad = ref)

  relabel <- c("causal", "reactive", "independent", "undecided")
  names(relabel) <- c("m.d_t.m", "t.d_m.t", "t.d_m.d", "t.md_m.d")
  best$triad <- factor(relabel[best$triad], relabel)
  best$alt <- factor(relabel[best$alt], relabel)

  result <- dplyr::arrange(
    dplyr::rename(
      dplyr::left_join(best, annotation, by = "pheno"),
      symbol = pheno),
    pv)

  attr(result, "pos") <- pos_t
  attr(result, "lod") <- lod_t
  attr(result, "target") <- colnames(target)

  class(result) <- c("mediate1_test", class(result))
  result
}
