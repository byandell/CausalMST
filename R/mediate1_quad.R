# Mediation quatrad tests
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
#' @param test 
#' @param pos Position of driver.
#' @param ... additional parameters
#'
#' @importFrom purrr map transpose
#' @importFrom stringr str_replace
#' @importFrom qtl2scan fit1 get_common_ids
#' @importFrom ggplot2 aes autoplot facet_wrap geom_hline geom_point ggplot
#'
#' @export
#'
mediate1_quad <- function(driver, target, mediator, fitFunction,
                          kinship=NULL, cov_tar=NULL, cov_med=NULL,
                          annotation, test = c("wilc","binom","joint","norm"),
                          pos = NULL,
                          ...) {

  test <- match.arg(test)
  
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
  
  quatrad_scan <- function(x, driver) {
    # Force x (= mediator column) to be matrix.
    x <- as.matrix(x)
    rownames(x) <- rownames(driver)
    colnames(x) <- "M"
    # Fit mediation models.
    models_par <- mediationModels(driver, target, x, 
                                  qtl2scan::fit1,
                                  kinship, cov_tar, cov_med,
                                  common = TRUE)
    # CMST on quatrads
    out <- quatrad_CMST(models_par$models, test, threshold=0.1)
    # Mediation LOD
    med_lod <- sum(models_par$comps$LR[c("t.d_t", "mediation")]) / log(10)
    out$mediation <- med_lod
    
    out
  }

  best <- purrr::map(as.data.frame(mediator), 
                     quatrad_scan, 
                     driver)
  best <- dplyr::bind_rows(best,
                          .id = "id")
  
  tmpfn <- function(x, pv, best.pv) {
    paste(unique(x[pv == min(pv) & best.pv == min(best.pv[pv == min(pv)])]),
          collapse = ",")
  }
  best <- dplyr::ungroup(
    dplyr::summarize(
      dplyr::group_by(
        dplyr::filter(best, !is.na(best.pv)),
        id),
      role = tmpfn(role, pv, best.pv),
      ref = tmpfn(ref, pv, best.pv),
      pv = min(pv),
      mediation = mediation[1]))

  result <- dplyr::left_join(best, annotation, by = "id")
  node_id <- quatrads()$node_id
  result <- dplyr::mutate(result,
                          role = factor(role, names(node_id)))
  
  attr(result, "pos") <- pos_t
  attr(result, "lod") <- lod_t
  
  class(result) <- c("mediate1_test", class(result))
  result
}
