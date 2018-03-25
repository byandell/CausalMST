#' Fit target relative to driver
#'
#' Fit a model for a target and get detailed results
#' about estimated coefficients and individuals contributions to the LOD score.
#' 
#' @details If \code{kinship} is absent, regression is performed.
#' If \code{kinship} is provided, a linear mixed model is used, with a
#' random effect estimated under the null hypothesis of no driver,
#' and then taken as fixed and known with driver.
#' The default version ignores kinship. See \code{\link[qtl2]{fit1}}
#' for use of \code{kinship}.
#'
#' @param driver A matrix of drivers, individuals x drivers
#' @param target A numeric vector of target values
#' @param kinship Optional kinship matrix.
#' @param addcovar An optional matrix of additive covariates.
#' @param nullcovar An optional matrix of additional additive
#' covariates that are used under the null hypothesis (of no driver)
#' but not under the alternative (with driver).
#' @param intcovar An optional matrix of interactive covariates.
#' @param weights An optional vector of positive weights for the
#' individuals. As with the other inputs, it must have `names`
#' for individual identifiers. Ignored if `kinship` is provided.#' 
#' 
#' @return A list containing
#' * `LR` - The overall likelihood ratio.
#' * `indLR` - Vector of individual contributions to the likelihood ratio.
#' * `df` - Model degrees of freedom.
#' 
#' @export
fitDefault <- function(driver,
                  target,
                  kinship = NULL,
                  addcovar = NULL, nullcovar=NULL,
                  intcovar=NULL, weights=NULL) {
  ll1 <- GetLogLik(driver, target, addcovar, intcovar)
  ll0 <- GetLogLik(NULL, target, cbind(addcovar, nullcovar))
  
  list(LR = ll1$logLik - ll0$logLik,
       indLR = ll1$ind_logLik - ll0$ind_logLik,
       df = ll1$df - ll0$df)
}