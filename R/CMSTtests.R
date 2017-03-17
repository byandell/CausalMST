#' CMST Tests for cross object
#' 
#' Set up and run CMST tests for cross object. See \code{\link{cmst}}.
#' 
#' @param cross object of class cross; see \code{\link[qtl]{read.cross}}.
#' @param pheno1 first phenotype name
#' @param pheno2 other phenotype names to test with \code{pheno1}
#' @param Q.chr,Q.pos chromosome name (character) and position (in map units).
#' @param addcov1,addcov2,intcov1,intcov2 additive and interactive covariates for \code{pheno1} and \code{pheno2}.
#' @param method method for CMST test (parametric, non-parametric, joint or all three); can provide more than one value.
#' @param penalty type of information criteria penalty (for BIC or AIC)
#' @param verbose verbose output if \code{TRUE}
#' 
#' @export
#' 
CMSTtests <- function(cross, 
                      pheno1, 
                      pheno2,
                      Q.chr,
                      Q.pos,
                      addcov1 = NULL, 
                      addcov2 = NULL, 
                      intcov1 = NULL, 
                      intcov2 = NULL, 
                      method = c("all", "par", "non.par", "joint"),
                      penalty = c("both", "bic", "aic"),
                      verbose = FALSE)
{
  if(length(pheno2) > 1)
    return(CMSTtestsList(cross, pheno1, pheno2, Q.chr, Q.pos,
                         addcov1, addcov2, intcov1, intcov2, 
                         method, penalty, verbose))
  
  setup <- qtl2cmst(cross, pheno1, pheno2,
                    Q.chr, Q.pos,
                    addcov1, addcov2, 
                    intcov1, intcov2)
  
  cmst(setup$pheno, setup$geno, 
       setup$resp_names, 
       setup$addcov, setup$intcov, 
       method, penalty)
}
CMSTtestsList <- function(cross, 
                          pheno1, 
                          pheno2,
                          Q.chr,
                          Q.pos,
                          addcov1 = NULL, 
                          addcov2 = NULL, 
                          intcov1 = NULL, 
                          intcov2 = NULL, 
                          method = c("par", "non.par", "joint", "all"),
                          penalty = c("bic", "aic", "both"),
                          verbose = TRUE)
{
  if(length(pheno2) == 1)
    return(CMSTtests(cross, pheno1, pheno2, Q.chr, Q.pos,
                     addcov1, addcov2, intcov1, intcov2, 
                     method, penalty, verbose))
  
  setup <- qtl2cmst(cross, pheno1, pheno2,
                    Q.chr, Q.pos,
                    addcov1, addcov2, 
                    intcov1, intcov2)
  
  cmsts(setup$pheno, setup$geno, 
        setup$resp_names, 
        setup$addcov, setup$intcov, 
        method, penalty)
}


