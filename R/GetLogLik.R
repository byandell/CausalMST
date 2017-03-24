GetLogLik <- function(driver, outcome, addcov, intcov, ll_function = logLik_calcs, ...) {

  # Construct design matrix X
  if(is.null(driver))
    driver <- matrix(1, length(outcome), 1)
  
  # Form model matrix from additive covariates.
  if(!is.null(addcov)) {
    form <- as.formula(paste(" ~ ", paste(colnames(addcov), collapse = "+")))
    X <- model.matrix(form, data = addcov)[,-1]
    if(is.null(dim(X))) {
      X <- as.matrix(X)
    }
  } else 
    X <- NULL
  
  # Add interactive covaraties.
  if(!is.null(intcov)) {
    if(ncol(driver) > 1) {
      int.sub.matrix <- 
        model.matrix(as.formula(paste("~", colnames(intcov))), intcov)[,-1]
      driverbyintcov <- driver[,-1] * int.sub.matrix
      X <- cbind(driver, X, driverbyintcov)
    } else {
      X <- cbind(driver, X)
    }
  } else {
    if(!is.null(addcov)){
      X <- cbind(driver, X)
    } else {
      X <- driver
    }
  }
  
  # Calculate log likelihood components
  ll_function(outcome, X, ...)
}

