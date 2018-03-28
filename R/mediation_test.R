# Mediation tests
#
#' Develop mediation models from driver, target and mediator
#'
#' @param target vector or 1-column matrix with target values
#' @param mediator matrix of mediators
#' @param annotation optional annotation data frame for mediators
#' @param driver vector or matrix with driver values
#' @param cov_tar optional covariates for target
#' @param cov_med optional covariates for mediator
#' @param kinship optional kinship matrix among individuals
#' @param driver_med optional driver matrix for mediators
#' @param test Type of CMST test.
#' @param pos Position of driver.
#' @param fitFunction function to fit models with driver, target and mediator
#' @param data_type Type of mediator data.
#' @param ... additional parameters
#'
#' @importFrom purrr map transpose
#' @importFrom stringr str_replace
#' @importFrom qtl2 decomp_kinship fit1 get_common_ids
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
mediation_test <- function(target, mediator, annotation, driver,
                          cov_tar=NULL, cov_med=NULL, kinship=NULL,
                          driver_med = NULL,
                          test = c("wilcoxon","binomial","joint","normal"),
                          pos = NULL,
                          fitFunction = fitQtl2,
                          data_type = c("phenotype","expression"),
                          ...) {
  
  ## Need to enable different covariates for different mediators.

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
                   wilcoxon = wilcIUCMST,
                   binomial = binomIUCMST,
                   joint    = normJointIUCMST,
                   normal   = normIUCMST)

  pos_tar <- pos
  
  # Make sure covariates are numeric
  cov_tar <- covar_df_mx(cov_tar)

  scan_max <- fitFunction(driver, target, kinship, cov_tar)
  LR_tar <- scan_max$LR
  
  use_1_driver <- is.null(annotation$driver) | is.null(driver_med)
  if(use_1_driver & !is.null(driver_med))
    driver_med <- NULL
  
  commons <- common_data(target, mediator, driver,
                         cov_tar, NULL, kinship,
                         common = use_1_driver)
  if(is.null(commons))
    return(NULL)
  
  driver <- commons$driver
  target <- commons$target
  kinship <- commons$kinship
  cov_tar <- commons$cov_tar
  common <- commons$common

  # Two reasons not to put cov_med in common_data call:
  # 1: different mediators may have different covariates
  # 2: cov_med is data frame, so need to be careful.
  # Fix up cov_med to match the rest
  m <- match(rownames(driver), rownames(cov_med), nomatch = 0)
  cov_med <- cov_med[m,, drop = FALSE]
  
  # Reorganize annotation and mediator data.
  # Need to make sure elements of mediator have same ids.
  mediator <- as.data.frame(commons$mediator)
  annotation <- dplyr::filter(
    annotation,
    id %in% colnames(mediator))

  # Workhorse: CMST on each mediator.
  best <- purrr::map(
    purrr::transpose(list(mediator = mediator,
                          annotation = purrr::transpose(annotation))),
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
    params = list(pos = pos_tar,
                  LR = LR_tar,
                  target = colnames(target),
                  data_type = data_type),
    targetFit = scan_max)
    
  class(result) <- c("mediation_test", class(result))
  result
}
#' @export
subset.mediation_test <- function(object, not_type, ...) {
  attrc <- class(object)
  object$best <- dplyr::filter(object$best, 
                               biotype != not_type)
  class(object) <- attrc
  
  object
}
