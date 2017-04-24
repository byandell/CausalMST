common_data <- function(driver, target, mediator,
                        kinship=NULL, cov_tar=NULL, cov_med=NULL) {
  
  # Make sure all are matrices
  target <- as.matrix(target)
  if(is.null(colnames(target)))
    colnames(target) <- "T"
  if(is.null(rownames(target)))
    rownames(target) <- seq_len(nrow(target))
  
  driver <- as.matrix(driver)
  if(is.null(rownames(driver)))
    rownames(driver) <- rownames(target)
  if(is.null(colnames(driver)))
    colnames(driver) <- "D"
  
  mediator <- as.matrix(mediator)
  if(is.null(rownames(mediator)))
    rownames(mediator) <- rownames(target)
  if(is.null(colnames(mediator)))
    colnames(mediator) <- "M"
  
  if(!is.null(cov_tar)) {
    cov_tar <- as.matrix(cov_tar)
    if(is.null(colnames(cov_tar)))
      colnames(cov_tar) <- paste0("covT", seq_len(ncol(cov_tar)))
  }
  if(!is.null(cov_med)) {
    cov_med <- as.matrix(cov_med)
    if(is.null(colnames(cov_med)))
      colnames(cov_med) <- paste0("covM", seq_len(ncol(cov_med)))
  }
  
  # Keep individuals with full records.
  ind2keep <-
    qtl2scan::get_common_ids(driver, target, mediator, cov_tar, cov_med, kinship,
                             complete.cases = TRUE)
  driver <- driver[ind2keep,, drop = FALSE]
  target <- target[ind2keep,, drop = FALSE]
  mediator <- mediator[ind2keep,, drop = FALSE]
  if(!is.null(cov_tar))
    cov_tar <- cov_tar[ind2keep,, drop = FALSE]
  if(!is.null(cov_med))
    cov_med <- cov_med[ind2keep,, drop = FALSE]
  if(!is.null(kinship)) {
    kinship <- kinship[ind2keep, ind2keep]
    kinship <- qtl2scan::decomp_kinship(kinship)
  }
  list(driver = driver,
       target = target,
       mediator = mediator,
       kinship = kinship, 
       cov_tar = cov_tar, 
       cov_med = cov_med)
}
