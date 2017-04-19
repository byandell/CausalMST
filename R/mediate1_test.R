# Mediation tests
#
#' Develop mediation models from driver, target and mediator
#'
#' @param driver vector or matrix with driver values
#' @param target vector or 1-column matrix with target values
#' @param mediator matrix with mediator values by column
#' @param fitFunction function to fit models with driver, target and mediator
#' @param kinship optional kinship matrix among individuals
#' @param cov_tar optional covariates for target
#' @param cov_med optional covariates for mediator
#' @param annotation Table with annotation, with \code{id} 
#' agreeing with column names of \code{mediator}.
#' @param test 
#' @param pos Position of driver.
#' @param lod_threshold LOD threshold to include mediator.
#' @param ... additional parameters
#'
#' @importFrom purrr map transpose
#' @importFrom stringr str_replace
#' @importFrom qtl2scan fit1 get_common_ids
#' @importFrom ggplot2 aes autoplot facet_grid geom_hline geom_point ggplot
#'
#' @export
#'
mediate1_test <- function(driver, target, mediator, fitFunction,
                          kinship=NULL, cov_tar=NULL, cov_med=NULL,
                          annotation, test = c("wilc","binom","joint","norm"),
                          pos = NULL,
                          lod_threshold = 5.5,
                          ...) {

  test <- match.arg(test)
  testfn <- switch(test,
                   wilc = CausalMST::wilcIUCMST,
                   binom = CausalMST::binomIUCMST,
                   joint = CausalMST::normJointIUCMST,
                   norm = CausalMST::normIUCMST)
  tmpfn <- function(x, models) {
    models <- subset(models, x)
    dplyr::filter(
      testfn(models),
      pv == min(pv))
  }
  
  pos_t <- pos

  scan_max <- fitFunction(driver, target, kinship, cov_tar)
  lod_t <- scan_max$lod
  
  commons <- common_data(driver, target, mediator,
                         kinship, cov_tar, cov_med)
  driver <- commons$driver
  target <- commons$target
  mediator <- commons$mediator
  kinship <- commons$kinship
  cov_tar <- commons$cov_tar
  cov_med <- commons$cov_med
  
  cmst_fit <- function(x, driver) {
    # Force x (= mediator column) to be matrix.
    x <- as.matrix(x)
    rownames(x) <- rownames(driver)
    colnames(x) <- "M"
    # Fit mediation models.
    models_par <- mediationModels(driver, target, x, 
                                  qtl2scan::fit1,
                                  kinship, cov_tar, cov_med,
                                  common = TRUE)
    # CMST on quatrads
    out <- dplyr::filter(
      testfn(subset(models_par$models, 1:4)),
      pv == min(pv))
    # Mediation LOD
    med_lod <- sum(models_par$comps$LR[c("t.d_t", "mediation")]) / log(10)
    # Mediator LOD
    medor_lod <- models_par$comp$LR["m.d_m"] / log(10)
    out$mediation <- med_lod
    out$mediator <- medor_lod
    
    out
  }

  best <- purrr::map(as.data.frame(mediator), 
                     cmst_fit, 
                     driver)
  best <- dplyr::rename(
    dplyr::filter(
      dplyr::bind_rows(best, .id = "id"),
      mediator >= lod_threshold),
    triad = ref)
  
  relabel <- c("D->M->T", "D->T->M", "M<-D->T", "D->{M,T}")
  names(relabel) <- c("m.d_t.m", "t.d_m.t", "t.d_m.d", "t.md_m.d")
  best$triad <- relabel[best$triad]

  result <- dplyr::left_join(best, annotation, by = "id")
  node_id <- quatrads()$node_id

  attr(result, "pos") <- pos_t
  attr(result, "lod") <- lod_t
  attr(result, "target") <- colnames(target)
  
  class(result) <- c("mediate1_test", class(result))
  result
}

#' @export
plot.mediate1_test <- function(x, ...)
  ggplot2::autoplot(x, ...)
#' @export
autoplot.mediate1_test <- function(x, ...)
  plot_mediate1_test(x, ...)
#' @export
plot_mediate1_test <- function(x, type = c("pos_lod","pos_pv","pv_lod"),
                               main = attr(x, "target"), ...) {
  type <- match.arg(type)
  
  pos_t <- attr(x, "pos")
  lod_t <- attr(x, "lod")
  
  switch(type,
         pos_pv = {
           p <- ggplot2::ggplot(x, 
               ggplot2::aes(x=pos, y=-log10(pv), col=triad, symbol=symbol)) +
             geom_point() +
             facet_grid(~triad)
           if(!is.null(pos_t))
             p <- p +
               geom_vline(xintercept = pos_t, col = "darkgrey")
         },
         pv_lod = {
           p <- ggplot2::ggplot(x, 
             ggplot2::aes(y=mediation, x=-log10(pv), col=triad, symbol=symbol)) +
           geom_point() +
           facet_grid(~triad) +
           geom_hline(yintercept = lod_t, col = "darkgrey")
         },
         pos_lod = {
           p <- ggplot2::ggplot(x, 
               ggplot2::aes(y=mediation, x=pos, col=triad, symbol=symbol)) +
             geom_point() +
             geom_hline(yintercept = lod_t, col = "darkgrey")
           if(!is.null(pos_t))
             p <- p +
               geom_vline(xintercept = pos_t, col = "darkgrey")
         })
  p + ggtitle(main)
}
