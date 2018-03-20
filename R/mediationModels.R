# Mediation plan
#
# scatter plot of mediator vs target with colors for driver and/or covariates
# modify driver position at will
# mediation idea with module of co-mapping traits
# how to address complex traits (multiple QTLs) with mediation
#
#' Develop mediation models from driver, target and mediator
#'
#' @param driver vector or matrix with driver values
#' @param target vector or 1-column matrix with target values
#' @param mediator vector or 1-column matrix with mediator values
#' @param fitFunction function to fit models with driver, target and mediator
#' @param ... additional parameters, which might include
#' @param kinship optional kinship matrix among individuals
#' @param cov_tar optional covariates for target
#' @param cov_med optional covariates for mediator
#'
#' @importFrom purrr map transpose
#' @importFrom stringr str_replace
#' @importFrom qtl2 fit1 get_common_ids
#'
#' @export
#'
mediationModels <- function(driver, target, mediator,
                     fitFunction = fitFunction, ...) {

  fits <- med_fits(driver, target, mediator, fitFunction, ...)

  list(models = fit_models(fits),
       comps  = fit_comps(fits))
}
fit_models <- function(fits, combos = combo_models()) {

  models <-
    purrr::transpose(
      purrr::map(combos,
                 combine_models, fits[c("lod", "ind_lod", "df")]))

  # Change LOD to LR
  names(models) <- stringr::str_replace(names(models), "lod", "LR")
  names(models) <- stringr::str_replace(names(models), "_LR", "LR")
  models$LR <- unlist(models$LR) * log(10)
  models$indLR <- as.data.frame(models$indLR) * log(10)
  
  models$df <- unlist(models$df)
  models$coef <- fits$coef

  class(models) <- c("cmst_models", class(models))

  models
}
combo_models <- function() {
  combos <- matrix(0, 5, 4,
                   dimnames = list(c("t.d_t", "t.md_t.m", "m.d_m", "t.m_t", "m.t_m"),
                                   c("t.d_m.t", "m.d_t.m", "t.d_m.d", "t.md_m.d")))
  combos[c(1,5), 1] <- 1 # causal
  combos[c(3,4), 2] <- 1 # reactive
  combos[c(1,3), 3] <- 1 # independent
  combos[   2:4, 4] <- 1 # correlated
  as.data.frame(combos)
}
combine_models <- function(combos, fits) {
  
  list(lod = sum(fits$lod * combos),
       ind_lod = fits$ind_lod %*% combos,
       df = sum(fits$df * combos),
       coef = fits$coef)
}

fit_comps <- function(fits, combos = combo_comps()) {
  
  comps <- list(LR = apply(fits$lod * combos, 2, sum) * log(10))

  class(comps) <- c("cmst_models", class(comps))
  
  comps
}
combo_comps <- function() {
  combos <- matrix(0, 5, 3,
                   dimnames = list(c("t.d_t", "t.md_t.m", "m.d_m", "t.m_t", "m.t_m"),
                                   c("t.d_t", "m.d_m", "mediation")))
  combos[  1, 1] <- 1        # target contrast
  combos[  3, 2] <- 1        # mediator contrast
  combos[1:2, 3] <- c(-1, 1) # mediation contrast
  as.data.frame(combos)
}
