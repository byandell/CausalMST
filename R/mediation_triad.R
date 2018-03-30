#' Scatter plot for mediator and target
#' 
#' @param target vector or 1-column matrix with target values
#' @param mediator vector or 1-column matrix with mediator values
#' @param driver vector or matrix with driver values
#' @param covar_tar optional covariates for target
#' @param covar_med optional covariates for mediator
#' @param kinship optional kinship matrix among individuals
#' @param fitFunction function to fit models with driver, target and mediator
#' @param sdp SNP distribution pattern for plot colors
#' @param allele Driver has alleles if \code{TRUE}, otherwise allele pairs.
#' 
#' @export
#' 
#' @importFrom stringr str_split
#' @importFrom ggplot2 aes autoplot facet_wrap geom_hline geom_smooth 
#' geom_text ggplot ggtitle scale_color_discrete xlab ylab
#' 
mediation_triad <- function(target, mediator, driver,
                        covar_tar, covar_med, kinship, 
                        fitFunction,
                        sdp,
                        allele = TRUE) {
  
  # Make sure covariates are numeric
  covar_tar <- covar_df_mx(covar_tar)
  covar_med <- covar_df_mx(covar_med)

  commons <- common_data(target, mediator, driver, 
                         covar_tar, covar_med, kinship)
  
  cov_names <- colnames(covar_med)[!(colnames(covar_med) %in% colnames(covar_tar))]
  
  for(i in c("target","mediator"))
    colnames(commons[[i]]) <- i

  dat <- data.frame(commons$driver, commons$target, commons$mediator,
                    commons$covar_tar, commons$covar_med[,cov_names, drop = FALSE])
  
  genos <- colnames(commons$driver)
  if(allele)
    haplos <- genos
  else
    haplos <- unique(unlist(stringr::str_split(genos, "")))
  
  # Would like to have option to have line per haplo.
  # But that requires some regression style approach, such as dividing up data
  # or fitting allele model. Guess is this would involve fitting allele model,
  # getting estimates of slopes for each ellele interacted with mediator,
  # creating data frame, and adding this to ggplot object.
  # Probably signal this with sdp = NULL option?
  
  alt <- haplos[sdp_to_logical(sdp, haplos)]
  if(allele) {
    dat$geno <- apply(round(2 * commons$driver), 1, 
                      function(x)
                        paste(rep(genos, x), collapse = ""))
    dat$alt <- factor(round(2 * apply(dat[, alt, drop = FALSE], 1, sum)))
  } else {
    dat$geno <- apply(commons$driver, 1, 
                      function(x) genos[which.max(x)])
    dat$alt <- sapply(stringr::str_split(genos, ""), 
                      function(x, alt) sum(x %in% alt), alt)
  }
  dat$Sex <- c("Female", "Male")[1 + dat$sex]
  
  # Fit target and target|mediator models
  fit <- med_fits(driver, target, mediator,
                  fitFunction, kinship, covar_tar, covar_med)
  for(i in names(fit$coef)[1:2]) {
    tmp <- fit$coef[[i]][1:8]
    dat[[i]] <- c(as.matrix(dat[names(tmp)]) %*% tmp)
  }
  
  class(dat) <- c("mediation_triad", class(dat))
  
  dat
}
#' @param x object of class \code{mediation_triad}
#' @param \dots additional parameters for plotting
#' @param type type of plot: one of \code{("by_mediator", "by_target", "driver_offset", "driver")}
#' @param tname target name (default \code{"target"})
#' @param mname mediator name (default \code{"mediator"})
#' @param dname driver name (default \code{"driver"})
#' @param centerline horizontal line at value (default = \code{0}); set to \code{NA} for no line or \code{NULL} for mean
#' @param main main title (defautl \code{tname})
#' 
#' @rdname mediation_triad
#' @export
ggplot_mediation_triad <- function(x, ..., 
                             type = c("by_mediator", "by_target", "driver_offset", "driver"),
                             tname = "target", mname = "mediator", dname = "driver",
                             centerline = 0,
                             main = tname) {
  
  type <- match.arg(type)
  
  p <- ggplot2::ggplot(x) +
    ggplot2::aes(label = geno, col = alt) +
    ggplot2::scale_color_discrete(name = dname)

  switch(type,
         by_mediator = {
           p <- p + 
             ggplot2::aes(mediator, target) +
             ggplot2::xlab(mname) +
             ggplot2::ylab(tname)
         },
         by_target = {
           p <- p + 
             ggplot2::aes(target, mediator) +
             ggplot2::xlab(tname) +
             ggplot2::ylab(mname)
         },
         driver_offset = {
           p <- p + 
             ggplot2::aes(mediator, t.d_t - t.md_t.m) +
             ggplot2::xlab(mname) +
             ggplot2::ylab(paste(dname, "effect offset"))
         },
         driver = {
           p <- p + 
             ggplot2::aes(t.d_t, t.d_t - t.md_t.m) +
             ggplot2::xlab(paste(dname, "effect")) +
             ggplot2::ylab(paste(dname, "effect offset"))
         })
  if(is.null(centerline)) {
    switch(type,
           by_mediator = {
             centerline <- mean(x$target, na.rm = TRUE)
           },
           by_target = {
             centerline <- mean(x$mediator, na.rm = TRUE)
           },
           driver_offset, driver = {
             centerline <- mean(x$t.d_t - x$t.md_t.m, na.rm = TRUE)
           })
  }
  
  if(!is.na(centerline)) {
    p <- p +
      ggplot2::geom_hline(yintercept = centerline)
  }
    
  p + 
    ggplot2::geom_smooth(method = "lm", se=FALSE) +
    ggplot2::geom_text(size=3) +
    ggplot2::facet_wrap(~ Sex) +
    ggplot2::ggtitle(main)
}
#' @export
autoplot.mediation_triad <- function(x, ...) {
  ggplot_mediation_triad(x, ...)
}

# from qtl2pattern
sdp_to_logical <- function(sdp, haplos = LETTERS[1:8]) {
  sapply(sdp, function(x, haplos) {
    as.logical(intToBits(x)[seq_along(haplos)])
  }, haplos)
}

