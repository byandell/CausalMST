CreateDesignMatrix <- function(geno, pheno, addcov.nms, intcov.nms) {
  if(is.null(geno))
    geno <- matrix(1, nrow(pheno), 1)
  
  # Form model matrix from covariates.
  uname <- unique(c(addcov.nms, intcov.nms))
  if(!is.null(uname)) {
    pheno <- pheno[, uname, drop = FALSE]
    form <- as.formula(paste(" ~ ", paste(uname, collapse = "+")))
    X <- model.matrix(form, data = pheno)[,-1]
    if(is.null(dim(X))) {
      X <- as.matrix(X)
    }
  } else 
    X <- NULL

  # Combine genotype and covaraties.
  if(!is.null(intcov.nms)) {
    if(ncol(geno) > 1) {
      int.sub.matrix <- model.matrix(as.formula(paste("~", intcov.nms)), 
                                      pheno[, intcov.nms, drop = FALSE])[,-1]
      genobyintcov <- geno[,-1] * int.sub.matrix
      X <- cbind(geno, X, genobyintcov)
    } else {
      X <- cbind(geno, X)
    }
  } else {
    if(!is.null(addcov.nms)){
      X <- cbind(geno, X)
    } else {
      X <- geno
    }
  }
  X
}
