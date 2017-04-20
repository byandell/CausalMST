med_scatter <- function(driver, target, mediator, fitFunction,
                        kinship, cov_tar, cov_med) {

  
  commons <- common_data(driver, target, mediator,
                         kinship, cov_tar, cov_med)
  
  cov_names <- names(cov_med)[!(names(cov_med) %in% names(cov_tar))]

  data.frame(commons$driver, commons$target, commons$mediator,
               commons$cov_tar, commons$cov_med[,cov_names, drop = FALSE])
}
