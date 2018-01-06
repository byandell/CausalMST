# Mediation tests
#
#' Develop mediation models from driver, target and mediator
#'
#' @param driver vector or matrix with driver values
#' @param target vector or 1-column matrix with target values
#' @param mediator List with mediators and annotation.
#' @param kinship optional kinship matrix among individuals
#' @param cov_tar optional covariates for target
#' @param cov_med optional covariates for mediator
#' @param test Type of CMST test.
#' @param pos Position of driver.
#' @param fitFunction function to fit models with driver, target and mediator
#' @param data_type Type of mediator data.
#' @param ... additional parameters
#'
#' @importFrom purrr map transpose
#' @importFrom stringr str_replace
#' @importFrom qtl2scan fit1 get_common_ids
#' @importFrom dplyr arrange bind_rows desc filter group_by left_join mutate one_of rename ungroup
#' @importFrom tidyr gather
#' @importFrom ggplot2 aes autoplot element_blank facet_grid facet_wrap 
#' geom_hline geom_point geom_vline ggplot 
#' ggtitle scale_color_manual scale_shape_manual theme xlab ylab
#' @importFrom grid grid.newpage pushViewport viewport grid.layout
#' @importFrom RColorBrewer brewer.pal
#'
#' @export
#'
mediate1_test <- function(mediator, driver, target,
                          kinship=NULL, cov_tar=NULL, cov_med=NULL, 
                          driver_med = NULL,
                          test = c("wilc","binom","joint","norm"),
                          pos = NULL,
                          fitFunction = qtl2scan::fit1,
                          data_type = c("expression","phenotype"),
                          ...) {

  if(is.null(mediator))
    return(NULL)
  
  data_type = match.arg(data_type)
  # Need following in annotation for each data_type
  #   id = identifier of mediator
  #   biotype = type of biological measurement
  cmstfn <- switch(data_type,
                   expression = cmst_default,
                   phenotype = cmst_pheno)
  
  test <- match.arg(test)
  testfn <- switch(test,
                   wilc = CausalMST::wilcIUCMST,
                   binom = CausalMST::binomIUCMST,
                   joint = CausalMST::normJointIUCMST,
                   norm = CausalMST::normIUCMST)

  pos_t <- pos
  
  scan_max <- fitFunction(driver, target, kinship, cov_tar)
  lod_t <- scan_max$lod
  
  use_1_driver <- is.null(mediator[[2]]$driver) | is.null(driver_med)
  if(use_1_driver & !is.null(driver_med))
    driver_med <- NULL
  
  commons <- common_data(driver, target, mediator[[1]],
                         kinship, cov_tar, cov_med,
                         common = use_1_driver)
  if(is.null(commons))
    return(NULL)
  
  driver <- commons$driver
  target <- commons$target
  kinship <- commons$kinship
  cov_tar <- commons$cov_tar
  cov_med <- commons$cov_med
  common <- commons$common
  
  # Reorganize annotation and mediator data.
  # Need to make sure elements of mediator have same ids.
  mediator[[1]] <- as.data.frame(commons$mediator)
  mediator[[2]] <- dplyr::filter(mediator[[2]],
                                 id %in% colnames(mediator[[1]]))
  annotation <- mediator[[2]]
  mediator[[2]] <- purrr::transpose(annotation)

  # Workhorse: CMST on each mediator.
  best <- purrr::map(
    purrr::transpose(mediator[1:2]),
    cmstfn, driver, target, 
    kinship, cov_tar, cov_med,
    driver_med,
    fitFunction, testfn, common)

  best <- dplyr::rename(
    dplyr::bind_rows(best, .id = "id"),
    triad = ref)
  
  # Kludge until I figure out why last level.
  relabel <- c("causal", "reactive", "independent", "correlated")
  names(relabel) <- c("m.d_t.m", "t.d_m.t", "t.d_m.d", "t.md_m.d")
  best$triad <- factor(relabel[best$triad], relabel)
  best$alt <- factor(relabel[best$alt], relabel)
  
  result <- list(
    best = dplyr::arrange(
      dplyr::left_join(best, annotation, by = "id"),
      pvalue),
    params = list(pos = pos_t,
                  lod = lod_t,
                  target = colnames(target),
                  data_type = data_type),
    targetFit = scan_max)
    
  class(result) <- c("mediate1_test", class(result))
  result
}
#' @export
subset.mediate1_test <- function(object, not_type, ...) {
  attrc <- class(object)
  object$best <- dplyr::filter(object$best, 
                               biotype != not_type)
  class(object) <- attrc
  
  object
}

cmst_default <- function(object, driver, target, 
                      kinship, cov_tar, cov_med,
                      driver_med,
                      fitFunction, testFunction,
                      common = TRUE) {
  
  # Force x (= mediator column) to be matrix.
  mediator <- as.matrix(object[[1]])
  rownames(mediator) <- rownames(driver)
  colnames(mediator) <- "mediator"
  if(!is.null(driver_med))
    driver_med <- driver_med[,, object[[2]]$driver]

  # Fit mediation models.
  models_par <- mediationModels(driver, target, mediator, 
                                fitFunction,
                                kinship, cov_tar, cov_med,
                                driver_med,
                                common = common)

  # CMST on key models.
  out <- dplyr::rename(
    dplyr::filter(
      testFunction(subset(models_par$models, 1:4)),
      pv == min(pv)),
    pvalue = pv)
  
  # Mediation LOD
  out$mediation <- sum(models_par$comps$LR[c("t.d_t", "mediation")]) / log(10)
  
  # Mediator LOD
  out$lod_med <- models_par$comp$LR["m.d_m"] / log(10)
  
  # Coefficients
  coef_target <- as.data.frame(t(models_par$models$coef$t.md_t.m[seq_len(ncol(driver))]))
  coef_mediator <- as.data.frame(t(models_par$models$coef$m.d_m[seq_len(ncol(driver))]))
  names(coef_mediator) <- paste0(names(coef_mediator), "_m")
  dplyr::bind_cols(out, coef_target, coef_mediator)
}

cmst_pheno <- function(object, driver, target, 
                       kinship, cov_tar, cov_med,
                       driver_med,
                       fitFunction, testFunction,
                       common = TRUE) {

  # Set up covariates 
  covi <- unlist(object[[2]][colnames(cov_med)])
  cov_medi <- cov_med[, covi, drop = FALSE]
  
  cmst_default(object, driver, target, 
               kinship, cov_tar, cov_medi,
               driver_med,
               fitFunction, testFunction,
               common)
}  



