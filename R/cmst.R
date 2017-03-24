#' Causal Model Selection Tests for outcomes given driver and covariates
#' 
#' Run CMST tests the causal direction between the first outcome \code{(y1)} and all other outcomes \code{(y2, ...)},
#' adjusting for the driver and covariates. The models tested (ignoring covariates and noise) are:
#' \itemize{
#' \item{\code{M1.y2.y1.z}} : \code{y2 <- y1 <- z}
#' \item{\code{M2.y1.y2.z}} : \code{y1 <- y2 <- z}
#' \item{\code{M3.y1.z.y2}} : \code{y1 <- z -> y2}
#' \item{\code{M4.y12.z}} : \code{y1 <- z -> y2} and \code{y1 <-> y2}
#' }
#' in which the arrows indicate direction of causality. For \code{M4}, the directionality between \code{y1} and \code{y2} is ambiguous.
#' The noise on outcomes \code{yi} is assumed to be normal.
#' The driver \code{z} is assumed to be causal for one or more outcomes \code{(y1, y2, ...)}, and may be categorical or continuous.
#' The covariates \code{X}, qualitative or quantitative, may act additively or be interactive with the driver.
#' The \code{driver} is the driver for outcomes, with covariates acting on outcomess.
#' The \code{resp_names}, \code{addcov} and \code{intcov} names must all be valid names in the \code{outcomes} data frame.
#' 
#' @param outcomes data frame with outcomes
#' @param driver data frame with drivers
#' @param resp_names response names
#' @param addcov,intcov additive and interactive covariate names for responses
#' @param method method for CMST test (parametric, non-parametric, joint or all three); can provide more than one value.
#' @param penalty type of information criteria penalty (for BIC or AIC)
#' @param verbose verbose output if \code{TRUE}
#' @param ll_function log likelihood calculation function; see \code{\link{logLik_calcs}}
#' @param ... possible additional arguments to \code{ll_function}
#' 
#' @export
#' @importFrom dplyr tbl_df
#' @importFrom mnormt pmnorm
#' @importFrom corpcor is.positive.definite make.positive.definite
#' 
cmst <- function(driver, outcomes, covariates=NULL, addcov=NULL, intcov=NULL,
                 method = c("par", "non.par", "joint", "all"),
                 penalty = c("bic", "aic", "both"),
                 verbose = FALSE,
                 ll_function = logLik_calcs,
                 ...) {
  
  if(any(is.na(c(driver))) | 
     any(is.na(c(outcomes))))
    stop("missing values not allowed")
  if(!is.null(covariates)) {
    if(any(is.na(c(covariates))))
      stop("missing values not allowed")
  }
  
  resp_names <- names(outcomes)
  if(length(resp_names) < 2)
    stop("need at least 2 response names")
  n_ind <- nrow(outcomes)
  n_drv <- nrow(driver)
  n_cov <- nrow(covariates)
  if(n_ind != n_drv)
    stop("driver, outcomes and covariates must have same number of rows")
  if(!is.null(covariates)) {
    if(n_ind != n_cov | n_drv != n_cov)
      stop("driver, outcomes and covariates must have same number of rows")

    if(!is.list(addcov)) {
      addcov <- list(addcov)
      if(length(addcov) == 1)
        addcov[[2]] <- NULL
    }
    if(!is.list(intcov)) {
      intcov <- list(intcov)
      if(length(intcov) == 1)
        intcov[[2]] <- NULL
    }
    for(i in 1:2) {
      if(!is.null(addcov[[i]])) {
        # Make sure addcov has all intcov.
        uname <- unique(c(addcov[[i]], intcov[[i]]))
        if(!all(uname %in% names(covariates)))
          stop("some of covariate names not in covariates data frame")
        addcov[[i]] <- covariates[, uname, drop = FALSE]
      }
      if(!is.null(intcov[[i]]))
        intcov[[i]] <- covariates[, intcov[[i]], drop = FALSE]
    }
  }

  method <- match.arg(method, several.ok = TRUE)
  if("all" %in% method)
    method <- c("par", "non.par", "joint")
  penalty <- match.arg(penalty, several.ok = TRUE)
  if("both" %in% penalty)
    penalty <- c("bic", "aic")
  
  ntests <- length(resp_names) - 1
  
  aux <- vector(mode = "list", length = ntests)
  nms <- names(aux) <- paste(resp_names[1], resp_names[-1], sep = "_")
  for(k in 1 : ntests) {
    if(verbose)
      cat("outcomes2 = ", k, "\n")   
    aux[[k]] <- cmst1(driver, outcomes[, c(1, 1 + k)], addcov, intcov, 
                      method, penalty, ll_function, ...)
  }
  
  models <- model_setup()
  
  # Rearrange results.

  out <- vector(mode = "list", length = 6)
  names(out) <- c("pvals", "AIC", "BIC", "AIC.Z", "BIC.Z", "R2")

  # Now populate the list
  out$pvals <- vector(mode = "list", length = length(nms))
  names(out$pvals) <- nms
  pval_names <- 
    paste(rep(c("BIC","AIC"), 3),
          rep(c("par","nonpar","joint"), each = 2),
          sep = ".")
  names(pval_names) <- 
    paste(rep("pvals", 6),
          rep(c("p","np","j"), each = 2),
          rep(c("BIC","AIC"), 3),
          sep = ".")
  out$pvals <- lapply(aux, function(x, pval_names) {
    pvals <- grep("pvals", names(x))
    out <- t(as.data.frame(x[pvals]))
    rownames(out) <- pval_names[rownames(out)]
    out
  }, pval_names)
  for(stat in c("AIC","BIC","AIC.Z","BIC.Z")) {
    out[[stat]] <- t(sapply(aux, function(x, stat) x[[stat]], stat))
    rownames(out[[stat]]) <- nms
  }
  out$R2 <- t(sapply(aux, function(x) x$R2))
  rownames(out$R2) <- nms
  
  class(out) <- c("cmst", class(out))
  out
}
#' @export
print.cmst <- function(x, ...) {
  cat("\nCausal Model Select Test p-values by outcome pair\n")
  for(outcome_pair in names(x$pvals)) {
    cat(outcome_pair, "\n")
    print(x$pvals[[outcome_pair]])
  }
  for(IC in c("AIC","BIC")) {
    if(!is.null(x$AIC)) {
      cat("\n", IC, "statistics\n")
      print(x[[IC]])
      print(x[[paste0(IC, ".Z")]])
    }
  }
  cat("\nExplained variation (R2)")
  print(x$R2)
  invisible()
}

cmst1 <- function(driver, outcomes, addcov=NULL, intcov=NULL,
                 method = c("par", "non.par", "joint", "all"),
                 penalty = c("bic", "aic", "both"), 
                 ll_function, ...) {
  
  resp_names <- names(outcomes)
  if(length(resp_names) != 2)
    stop("only two outcomes for cmst1 call")
  
  outcomes <- outcomes[resp_names]

  # log.lik.1 #
  dfcol <- function(x, y, i) as.data.frame(cbind(x, as.matrix(y[,i, drop = FALSE])))
  ll_list <- list()
  ll_list[["y1.z"]] <-   
    GetLogLik(driver, outcomes[[1]], 
              addcov[[1]], intcov[[1]], 
              ll_function, ...)
  ll_list[["y2.z"]] <-   
    GetLogLik(driver, outcomes[[2]], 
              addcov[[2]], intcov[[2]], 
              ll_function, ...)
  ll_list[["y1.y2"]] <-  
    GetLogLik(NULL,   outcomes[[1]], 
              dfcol(addcov[[1]], outcomes, 2), intcov[[1]], 
              ll_function, ...)
  ll_list[["y2.y1"]] <-  
    GetLogLik(NULL,   outcomes[[2]], 
              dfcol(addcov[[2]], outcomes, 1), intcov[[2]], 
              ll_function, ...)
  ll_list[["y2.y1z"]] <- 
    GetLogLik(driver, outcomes[[2]], 
              dfcol(addcov[[2]], outcomes, 1), intcov[[2]], 
              ll_function, ...)
  
  models <- model_setup()
  
  R2 <- rep(0,2)
  names(R2) <- models$first[1:2]
  for(i in 1:2) {
    TSS <- sum((outcomes[[resp_names[i]]] - 
                  mean(outcomes[[resp_names[i]]]))^2)
    R2[models$first[i]] <- 1 - 
      (ll_list[[models$first[i]]]$RSS / TSS)
  }
  
  loglik <- model.dim <- rep(0,4)
  vec.logLik <- matrix(0, n_ind, 4)
  names(loglik) <- names(model.dim) <- colnames(vec.logLik) <- models$model
  for(i in 1:4) {
    loglik[models$model[i]] <- 
      ll_list[[models$first[i]]]$log.lik + 
      ll_list[[models$second[i]]]$log.lik
    model.dim[models$model[i]] <- 2 + 
      ll_list[[models$first[i]]]$d + 
      ll_list[[models$second[i]]]$d
    vec.logLik[, models$model[i]] <- 
      ll_list[[models$first[i]]]$vec.log.lik + 
      ll_list[[models$second[i]]]$vec.log.lik
  }
  
  out <- list(resp = resp_names,
              addcov = addcov,
              intcov = intcov,
              n.ind = n_ind,
              loglik = loglik,  
              model.dim = model.dim, 
              R2 = R2,
              S.hat = calcShat(vec.logLik))
  
  # Set up information criteria.
  tmpfn <- function(x) {
    x <- t(x)
    x <- x[!is.na(x)]
    names(x) <- c("z.12", "z.13", "z.14", "z.23", "z.24", "z.34")
    x
  }
  if("bic" %in% penalty & ("par" %in% method | "joint" %in% method)) {
    out$BIC <- calcBICs(out$n.ind, model.dim, loglik)
    BIC.Z <- calcZ(out$S.hat, out$BIC, out$n.ind)
    out$BIC.Z <- tmpfn(BIC.Z)
  }
  if("aic" %in% penalty & ("par" %in% method | "joint" %in% method)) {
    out$AIC <- calcAICs(out$n.ind, model.dim, loglik)
    AIC.Z <- calcZ(out$S.hat, out$AIC, out$n.ind)
    out$AIC.Z <- tmpfn(AIC.Z)
  }
  
  if("par" %in% method) {
    if("bic" %in% penalty)
      out$pvals.p.BIC <- ParametricIUCMST(BIC.Z)
    if("aic" %in% penalty)
      out$pvals.p.AIC <- ParametricIUCMST(AIC.Z)
  }  
  if("non.par" %in% method) {
    if("bic" %in% penalty)
      out$pvals.np.BIC <- NonparametricIUCMST("bic", out$n.ind, model.dim, vec.logLik)
    if("aic" %in% penalty)
      out$pvals.np.AIC <- NonparametricIUCMST("aic", out$n.ind, model.dim, vec.logLik)
  }  
  if("joint" %in% method) {
    Cor.hat <- corHat(out$S.hat)
    if("bic" %in% penalty)
      out$pvals.j.BIC <- ParametricJointCMST(BIC.Z, Cor.hat)
    if("aic" %in% penalty)
      out$pvals.j.AIC <- ParametricJointCMST(AIC.Z, Cor.hat) 
  }
  pvals <- grep("pvals", names(out))
  out[pvals] <- lapply(out[pvals], function(x,m) {
    names(x) <- m
    x
  }, models$model)
  
  out
}
model_setup <- function() {
  models <- dplyr::tbl_df(matrix(c(
    "M1.y2.y1.z", "y1.z", "y2.y1",
    "M2.y1.y2.z", "y2.z", "y1.y2",
    "M3.y1.z.y2", "y2.z", "y1.z",
    "M4.y12.z",   "y1.z", "y2.y1z"),
    4, 3, byrow = TRUE))
  names(models) <- c("model","first","second")
  models
}