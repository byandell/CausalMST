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
  
  # Genotype matrix.
  geno <- qtlMatrix(cross, Q.chr, Q.pos)
  
  # Phenotype and covariate matrix.  
  pheno_names <- unique(unlist(c(resp_names, addcov, intcov)))
  pheno <- cross$pheno[, pheno_names]
  # Reduce to individuals with no missing data.
  to.drop <- as.vector(attr(na.omit(cross$pheno[, pheno_names]), 
                            "na.action"))
  if (!is.null(to.drop)) {
    cross <- subset(cross, ind = -to.drop)
  }
  list(pheno = pheno, 
       geno = geno, 
       resp_names = resp_names, 
       addcov = addcov, 
       intcov = intcov)
}