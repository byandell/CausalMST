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
#' @importFrom qtl2scan fit1 get_common_ids
#'
#' @export
#'
mediationModels <- function(driver, target, mediator,
                     fitFunction = fitFunction, ...) {

  fits <- med_fits(driver, target, mediator, fitFunction, ...)

  list(models = fit_models(fits),
       comps  = fit_models(fits, combo_comps()))
}
fit_models <- function(fits, combos = combo_models()) {

  models <-
    purrr::transpose(
      purrr::map(combos,
                 comb_models, fits[c("lod", "ind_lod", "df")]))

  names(models) <- stringr::str_replace(names(models), "lod", "LR")
  names(models) <- stringr::str_replace(names(models), "_LR", "LR")
  models$LR <- unlist(models$LR) * log(10)
  models$indLR <- as.data.frame(models$indLR) * log(10)
  models$df <- unlist(models$df)
  models$coef <- fits$coef

  class(models) <- c("cmst_models", class(models))

  models
}

med_fits <- function(driver, target, mediator, fitFunction,
                     kinship=NULL, cov_tar=NULL, cov_med=NULL,
                     driver_med = NULL,
                     common = FALSE, ...) {

  if(!common) {
    commons <- common_data(driver, target, mediator,
                           kinship, cov_tar, cov_med, driver_med)
    driver <- commons$driver
    target <- commons$target
    mediator <- commons$mediator
    kinship <- commons$kinship
    cov_tar <- commons$cov_tar
    cov_med <- commons$cov_med
    driver_med <- commons$driver_med
  }
  if(is.null(driver_med))
    driver_med <- driver
  
  fits <- list(
    t.d_t =
      fitFunction(driver,
           target,
           kinship,
           cov_tar),
    t.md_t.m =
      fitFunction(driver,
           target,
           kinship,
           cbind(cov_tar, mediator)),
    m.d_m =
      fitFunction(driver_med,
           mediator,
           kinship,
           cov_med),
    t.m_t =
      fitFunction(cbind(1, mediator),
           target,
           kinship,
           cov_tar),
    m.td_m.t =
      fitFunction(driver_med,
                  mediator,
                  kinship,
                  cbind(cov_med, target)),
    m.t_m =
      fitFunction(cbind(1, target),
                  mediator,
                  kinship,
                  cov_med))

  # Transpose list
  fits <- purrr::transpose(fits)

  fits$lod <- unlist(fits$lod)
  fits$ind_lod <- as.matrix(as.data.frame(fits$ind_lod))

  # Add model degrees of freedom
  ndt <- ncol(driver) - 1
  ndm <- ncol(driver_med) - 1
  nmed <- ncol(mediator)
  fits$df <- c(rep(ndt, 2), ndm, nmed, ndm, nmed)
  names(fits$df) <- names(fits$lod)

  fits
}

combo_models <- function() {
  combos <- matrix(0, 11, 6)
  colnames(combos) <- c("t.d_t", "t.md_t.m", "m.d_m", "t.m_t", "m.td_m.t", "m.t_m")
  combos[1, c(1,6)] <- 1
  combos[2, c(3,4)] <- 1
  combos[3, c(1,3)] <- 1
  combos[4, 2:4] <- 1
  combos[5, 1] <- 1
  combos[6, 3] <- 1
  combos[7, c(2,4)] <- 1
  combos[8, 5:6] <- 1
  combos[9, 4] <- 1
  combos[11, c(1,5:6)] <- 1
  combos <- as.data.frame(t(combos))
  names(combos) <-
    c("t.d_m.t", "m.d_t.m", "t.d_m.d", "t.md_m.d", "t.d_m",
      "m.d_t", "t.md_m", "m.td_t", "t.m_m", "t_m", "m.td_t.d")
  combos
}
combo_comps <- function() {
  combos <- matrix(0, 7, 6)
  colnames(combos) <- c("t.d_t", "t.md_t.m", "m.d_m", "t.m_t", "m.td_m.t", "m.t_m")
  combos[1, 1] <- 1
  combos[2, 2] <- 1
  combos[3, 3] <- 1
  combos[4, 5] <- 1
  combos[5, 4] <- 1
  combos[6, c(1:2,4)] <- c(-1, 1, 1)
  combos[7, 1:2] <- c(-1, 1)
  combos <- as.data.frame(t(combos))
  names(combos) <-
    c("t.d_t", "t.md_t.m", "m.d_m", "m.td_m.t", "t.m_t", "t.md_t.d", "mediation")
  combos
}
comb_models <- function(combos, fits) {

  list(lod = sum(fits$lod * combos),
       ind_lod = fits$ind_lod %*% combos,
       df = sum(fits$df * combos),
       coef = fits$coef)
}
