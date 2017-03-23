CreateDesignMatrix <- function(driver, covariates, addcov.nms, intcov.nms) {
  if(is.null(driver))
    driver <- matrix(1, nrow(covariates), 1)
  
  # Form model matrix from covariates.
  uname <- unique(c(addcov.nms, intcov.nms))
  if(!is.null(uname)) {
    covariates <- covariates[, uname, drop = FALSE]
    form <- as.formula(paste(" ~ ", paste(uname, collapse = "+")))
    X <- model.matrix(form, data = covariates)[,-1]
    if(is.null(dim(X))) {
      X <- as.matrix(X)
    }
  } else 
    X <- NULL

  # Combine driver and covaraties.
  if(!is.null(intcov.nms)) {
    if(ncol(driver) > 1) {
      int.sub.matrix <- model.matrix(as.formula(paste("~", intcov.nms)), 
                                      covariates[, intcov.nms, drop = FALSE])[,-1]
      driverbyintcov <- driver[,-1] * int.sub.matrix
      X <- cbind(driver, X, driverbyintcov)
    } else {
      X <- cbind(driver, X)
    }
  } else {
    if(!is.null(addcov.nms)){
      X <- cbind(driver, X)
    } else {
      X <- driver
    }
  }
  X
}
