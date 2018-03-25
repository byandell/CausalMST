cmst_default <- function(object, driver, target, 
                         kinship, cov_tar, cov_med,
                         driver_med,
                         fitFunction, testFunction,
                         common = TRUE) {
  
  # Force x (= mediator column) to be matrix.
  mediator <- as.matrix(object[[1]])
  rownames(mediator) <- rownames(driver)
  colnames(mediator) <- "mediator"
  if(!is.null(driver_med))
    driver_med <- driver_med[,, object[[2]]$driver]
  
  # Make sure covariates are numeric
  cov_med <- covar_df_mx(cov_med)
  
  # Fit mediation models.
  models_par <- mediationModels(driver, target, mediator, 
                                fitFunction,
                                kinship, cov_tar, cov_med,
                                driver_med,
                                common = common)
  
  # CMST on key models. Pick first as best.
  out <- head(dplyr::rename(
    dplyr::filter(
      testFunction(models_par$models),
      pv == min(pv)),
    pvalue = pv), n = 1L)
  
  # Mediation LR
  out$mediation <- sum(models_par$comps$LR[c("t.d_t", "mediation")])
  
  # Mediator LR
  out$LRmed <- models_par$comps$LR["m.d_m"]
  
  # Coefficients
  coef_target <- as.data.frame(t(models_par$models$coef$t.md_t.m[seq_len(ncol(driver))]))
  coef_mediator <- as.data.frame(t(models_par$models$coef$m.d_m[seq_len(ncol(driver))]))
  names(coef_mediator) <- paste0(names(coef_mediator), "_m")
  dplyr::bind_cols(out, coef_target, coef_mediator)
}

cmst_pheno <- function(object, driver, target, 
                       kinship, cov_tar, cov_med,
                       driver_med,
                       fitFunction, testFunction,
                       common = TRUE) {
  
  # Currently, mediate1_test uses elements of object[[2]] (columns of annot data frame)
  # to assess TRUE/FALSE on covariate columns. This will likely change.
  
  # Get covariate names appropriate for mediator 
  cov_names <- unlist(object[[2]][colnames(cov_med)])
  if(length(cov_names))
    cov_med <- cov_med[, cov_names, drop = FALSE]
  else
    cov_med <- NULL
  
  cmst_default(object, driver, target, 
               kinship, cov_tar, cov_med,
               driver_med,
               fitFunction, testFunction,
               common)
}  
