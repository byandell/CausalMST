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
                 comb_models, fits))

  names(models) <- stringr::str_replace(names(models), "lod", "LR")
  names(models) <- stringr::str_replace(names(models), "_LR", "LR")
  models$LR <- unlist(models$LR) * log(10)
  models$indLR <- as.data.frame(models$indLR) * log(10)
  models$df <- unlist(models$df)

  class(models) <- c("cmst_models", class(models))

  models
}

common_data <- function(driver, target, mediator,
                        kinship=NULL, cov_tar=NULL, cov_med=NULL) {

  # Make sure all are matrices
  target <- as.matrix(target)
  if(is.null(colnames(target)))
    colnames(target) <- "T"
  if(is.null(rownames(target)))
    rownames(target) <- seq_len(nrow(target))
  
  driver <- as.matrix(driver)
  if(is.null(rownames(driver)))
    rownames(driver) <- rownames(target)
  if(is.null(colnames(driver)))
    colnames(driver) <- "D"
  
  mediator <- as.matrix(mediator)
  if(is.null(rownames(mediator)))
    rownames(mediator) <- rownames(target)
  if(is.null(colnames(mediator)))
    colnames(mediator) <- "M"
  
  if(!is.null(cov_tar)) {
    cov_tar <- as.matrix(cov_tar)
    if(is.null(colnames(cov_tar)))
      colnames(cov_tar) <- paste0("covT", seq_len(ncol(cov_tar)))
  }
  if(!is.null(cov_med)) {
    cov_med <- as.matrix(cov_med)
    if(is.null(colnames(cov_med)))
      colnames(cov_med) <- paste0("covM", seq_len(ncol(cov_med)))
  }
  
  # Keep individuals with full records.
  ind2keep <-
    qtl2scan::get_common_ids(driver, target, mediator, cov_tar, cov_med, kinship,
                             complete.cases = TRUE)
  driver <- driver[ind2keep,, drop = FALSE]
  target <- target[ind2keep,, drop = FALSE]
  mediator <- mediator[ind2keep,, drop = FALSE]
  if(!is.null(cov_tar))
    cov_tar <- cov_tar[ind2keep,, drop = FALSE]
  if(!is.null(cov_med))
    cov_med <- cov_med[ind2keep,, drop = FALSE]
  if(!is.null(kinship)) {
    kinship <- kinship[ind2keep, ind2keep]
    kinship <- qtl2scan::decomp_kinship(kinship)
  }
  list(driver = driver,
       target = target,
       mediator = mediator,
       kinship = kinship, 
       cov_tar = cov_tar, 
       cov_med = cov_med)
}
med_fits <- function(driver, target, mediator, fitFunction,
                     kinship=NULL, cov_tar=NULL, cov_med=NULL,
                     common = FALSE, ...) {

  if(!common) {
    commons <- common_data(driver, target, mediator,
                           kinship, cov_tar, cov_med)
    driver <- commons$driver
    target <- commons$target
    mediator <- commons$mediator
    kinship <- commons$kinship
    cov_tar <- commons$cov_tar
    cov_med <- commons$cov_med
  }
  
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
      fitFunction(driver,
           mediator,
           kinship,
           cov_med),
    t.m_t =
      fitFunction(mediator,
           target,
           kinship,
           cov_tar))

  # Transpose list
  fits <- purrr::transpose(fits)

  fits$lod <- unlist(fits$lod)
  fits$ind_lod <- as.matrix(as.data.frame(fits$ind_lod))

  # Add model degrees of freedom
  fits$df <- c(rep(ncol(driver) - 1, 3), ncol(mediator))
  names(fits$df) <- names(fits$lod)

  fits
}

combo_models <- function() {
  combos <- matrix(0, 10, 4)
  colnames(combos) <- c("t.d_t", "t.md_t.m", "m.d_m", "t.m_t")
  combos[1, c(1,4)] <- 1
  combos[2, c(3,4)] <- 1
  combos[3, c(1,3)] <- 1
  combos[4, 2:4] <- 1
  combos[5, 1] <- 1
  combos[6, 3] <- 1
  combos[7, c(2,4)] <- 1
  combos[8, ] <- c(-1, rep(1, 3))
  combos[9, 4] <- 1
  combos <- as.data.frame(t(combos))
  names(combos) <-
    c("t.d_m.t", "m.d_t.m", "t.d_m.d", "t.md_m.d", "t.d_m",
      "m.d_t", "t.md_m", "m.td_t", "t.m_m", "t_m")
  combos
}
combo_comps <- function() {
  combos <- matrix(0, 7, 4)
  colnames(combos) <- c("t.d_t", "t.md_t.m", "m.d_m", "t.m_t")
  combos[1, 1] <- 1
  combos[2, 2] <- 1
  combos[3, 3] <- 1
  combos[4, 1:3] <- c(-1, 1, 1)
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
       df = sum(fits$df * combos))
}
