#' Set up data for cmst function
#' 
#' Set up R/qtl cross object and identifiers for use with \code{\link{cmst}}.
#' The \code{pheno} and \code{cov} arguments must be valid names of \code{cross$pheno}.
#' The \code{Q.chr} must be a valid chromosome identifier for \code{cross}.
#' 
#' @param cross object of class cross; see \code{\link[qtl]{read.cross}}.
#' @param pheno1 first phenotype name
#' @param pheno2 other phenotype names to test with \code{pheno1}
#' @param Q.chr,Q.pos chromosome name (character) and position (in map units).
#' @param addcov1,addcov2,intcov1,intcov2 additive and interactive covariates for \code{pheno1} and \code{pheno2}.
#' 
#' @export
#' 
qtl2cmst <- function(cross, 
                     pheno1, 
                     pheno2,
                     Q.chr,
                     Q.pos,
                     addcov1 = NULL, 
                     addcov2 = NULL, 
                     intcov1 = NULL, 
                     intcov2 = NULL) {
  
  if (!any(class(cross) == "cross")) 
    stop("Input should have class \"cross\".")
  
  resp_names <- c(pheno1, pheno2)
  addcov <- list(addcov1, addcov2)
  intcov <- list(intcov1, intcov2)
  unames <- unique(unlist(c(addcov, intcov)))
  
  # Reduce to individuals with no missing data.
  to.drop <- as.vector(attr(na.omit(cross$pheno[, c(resp_names,unames)]), 
                            "na.action"))
  if (!is.null(to.drop)) {
    cross <- subset(cross, ind = -to.drop)
  }

  # Genotype matrix.
  geno <- qtlMatrix(cross, Q.chr, Q.pos)
  # Check for missing genotypes, and drop individuals with missing values.
  to.drop <- apply(geno, 1, function(x) any(is.na(x)))
  if(any(to.drop)) {
    geno <- geno[!do.drop, ]
    cross <- subset(cross, ind = !to.drop)
  }
  
  # Outcomes and covariates matrices.  
  outcomes <- cross$pheno[, resp_names, drop = FALSE]
  if(!is.null(unames)) {
    covariates <- cross$pheno[, unames, drop = FALSE]
  } else {
    covariates <- NULL
  }
  
  list(driver = geno, 
       outcomes = outcomes, 
       covariates = covariates, 
       addcov = addcov, 
       intcov = intcov)
}