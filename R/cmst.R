#' CMST Tests for genotype and phenotypes with covariates
#' 
#' Run CMST tests with one QTL and 2 or more phenotypes.
#' The \code{geno} is the driver for phenotype responses, with covariates acting on responses.
#' The \code{resp_names}, \code{addcov} and \code{intcov} names must all be valid names in the \code{pheno} data frame.
#' 
#' @param pheno data frame with phenotypes and covariates
#' @param geno data frame with QTL genotypes
#' @param resp_names response names
#' @param addcov,intcov additive and interactive covariate names for responses
#' @param method method for CMST test (parametric, non-parametric, joint or all three); can provide more than one value.
#' @param penalty type of information criteria penalty (for BIC or AIC)
#' @param verbose verbose output if \code{TRUE}
#' 
#' @export
#' 
cmst <- function(pheno, geno, resp_names, addcov, intcov,
                 method = c("par", "non.par", "joint", "all"),
                 penalty = c("bic", "aic", "both"),
                 verbose = FALSE) {
  
  if(length(resp_names) < 2)
    stop("need at least 2 response names")
  if(nrow(pheno) != nrow(geno))
    stop("pheno and geno must have same number of rows")
  uname <- unique(c(resp_names, addcov, intcov))
  if(!all(uname %in% names(pheno)))
    stop("some of phenotype or covariate names not in pheno data frame")
  
  if(length(resp_names) > 2)
    return(cmsts(pheno, geno, resp_names, addcov, intcov,
                 method, penalty, verbose))
  
  method <- match.arg(method, several.ok = TRUE)
  if("all" %in% method)
    method <- c("par", "non.par", "joint")
  penalty <- match.arg(penalty, several.ok = TRUE)
  if("both" %in% penalty)
    penalty <- c("bic", "aic")
  
  # log.lik.1 #
  ll_list <- list()
  ll_list[["y1.q"]] <-   
    GetLogLik(geno, pheno, resp_names[1], addcov[[1]], intcov[[1]])
  ll_list[["y2.q"]] <-   
    GetLogLik(geno, pheno, resp_names[2], addcov[[2]], intcov[[2]])
  ll_list[["y1.y2"]] <-  
    GetLogLik(NULL, pheno, resp_names[1], c(addcov[[1]], resp_names[2]), intcov[[1]])
  ll_list[["y2.y1"]] <-  
    GetLogLik(NULL, pheno, resp_names[2], c(addcov[[2]], resp_names[1]), intcov[[2]])
  ll_list[["y2.y1g"]] <- 
    GetLogLik(geno, pheno, resp_names[2], c(addcov[[2]], resp_names[1]), intcov[[2]])
  
  models <- dplyr::tbl_df(matrix(c(
    "M1.y2.y1.q", "y1.q", "y2.y1",
    "M2.y1.y2.q", "y2.q", "y1.y2",
    "M3.y1.q.y2", "y2.q", "y1.q",
    "M4.y12.q",   "y1.q", "y2.y1g"),
    4, 3, byrow = TRUE))
  names(models) <- c("model","first","second")
  
  R2 <- rep(0,2)
  names(R2) <- models$first[1:2]
  for(i in 1:2) {
    TSS <- sum((pheno[[resp_names[i]]] - 
                  mean(pheno[[resp_names[i]]]))^2)
    R2[models$first[i]] <- 1 - 
      (ll_list[[models$first[i]]]$RSS / TSS)
  }

  loglik <- model.dim <- rep(0,4)
  vec.logLik <- matrix(0, nrow(pheno), 4)
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
              n.ind = nrow(pheno),
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
    out$BICs <- calcBICs(out$n.ind, model.dim, loglik)
    Z.bic <- calcZ(out$S.hat, out$BICs, out$n.ind)
    out$Z.bic <- tmpfn(Z.bic)
  }
  if("aic" %in% penalty & ("par" %in% method | "joint" %in% method)) {
    out$AICs <- calcAICs(out$n.ind, model.dim, loglik)
    Z.aic <- calcZ(out$S.hat, out$AICs, out$n.ind)
    out$Z.aic <- tmpfn(Z.aic)
  }
  
  if("par" %in% method) {
    if("bic" %in% penalty)
      out$pvals.p.BIC <- ParametricIUCMST(Z.bic)
    if("aic" %in% penalty)
      out$pvals.p.AIC <- ParametricIUCMST(Z.aic)
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
      out$pvals.j.BIC <- ParametricJointCMST(Z.bic, Cor.hat)
    if("aic" %in% penalty)
      out$pvals.j.AIC <- ParametricJointCMST(Z.aic, Cor.hat) 
  }
  out
}
cmsts <- function(pheno, geno, resp_names, addcov, intcov,
                  method, penalty, verbose = FALSE) {
  
  ntests <- length(resp_names) - 1
  
  aux <- vector(mode = "list", length = ntests)
  for(k in 1 : ntests) {
    if(verbose)
      cat("pheno2 = ", k, "\n")   
    aux[[k]] <- cmst(pheno, geno, 
                     resp_names[c(1, 1 + k)], 
                     addcov, intcov, 
                     method, penalty)
  }
  
  # Rearrange results.
  models <- dplyr::tbl_df(matrix(c(
    "M1.y2.y1.q", "y1.q", "y2.y1",
    "M2.y1.y2.q", "y2.q", "y1.y2",
    "M3.y1.q.y2", "y2.q", "y1.q",
    "M4.y12.q",   "y1.q", "y2.y1g"),
    4, 3, byrow = TRUE))
  names(models) <- c("model","first","second")
  
  nms <- paste(resp_names[1], resp_names[-1], sep = "_")
  pval.nms <- paste("pval", 1:4, sep = ".")
  z.names <- c("z.12", "z.13", "z.14", "z.23", "z.24", "z.34")
  AIC.nms <- c(paste("AIC", 1:4, sep = "."), z.names)
  BIC.nms <- c(paste("BIC", 1:4, sep = "."), z.names)
  
  out <- vector(mode = "list", length = 9)
  names(out) <- c("R2s", "AIC.stats", "BIC.stats", 
                  "pvals.j.BIC", "pvals.p.BIC", "pvals.np.BIC",
                  "pvals.j.AIC", "pvals.p.AIC", "pvals.np.AIC")

  # Now populate the list
  out$AIC.stats <- t(sapply(aux, function(x) c(x$AICs, x$Z.aic)))
  out$BIC.stats <- t(sapply(aux, function(x) c(x$BICs, x$Z.bic)))
  rownames(out$AIC.stats) <- rownames(out$BIC.stats) <- nms
  for(i in names(out)[c(1,4:9)]) {
    if(!is.null(aux[[1]][[i]])) {
      out[[i]] <- t(sapply(aux, function(x, i) x[[i]], i))
      dimnames(out[[i]]) <- list(nms, models$model)
    } else
      out[[i]] <- NULL
  }
  out
}