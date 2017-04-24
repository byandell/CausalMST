#' @export
plot.mediate1_test <- function(x, ...)
  ggplot2::autoplot(x, ...)
#' @export
autoplot.mediate1_test <- function(x, ...)
  plot_mediate1_test(x, ...)
#' @export
plot_mediate1_test <- function(x, type = c("pos_lod","pos_pvalue","pvalue_lod"),
                               main = params$target,
                               maxPvalue = 0.1, 
                               local_only = FALSE, 
                               significant = TRUE,
                               ...) {
  type <- match.arg(type)
  
  params <- attr(x, "params")
  pos_t <- params$pos
  lod_t <- params$lod
  main
  
  if(!("symbol" %in% names(x)))
    x <- dplyr::rename(x, symbol = id)
  
  relabel <- levels(x$triad)
  if(significant) {
    x <- dplyr::filter(x, x$pvalue <= maxPvalue)
  } else {
    relabel <- c(relabel, 
                 paste0("n.s. (p>", round(maxPvalue, 2), ")"))
    tmp <- as.character(x$triad)
    tmp[x$pvalue > maxPvalue] <- relabel[5]
    x$triad <- factor(tmp, levels = relabel)
  }
  x <- dplyr::arrange(x, dplyr::desc(triad))
  
  # For expression, use qtl_pos if not missing.
  if(params$data_type == "expression") {
    if(local_only)
      x <- dplyr::filter(x, local)
    else
      x <- dplyr::mutate(x, pos = ifelse(local, pos, qtl_pos))
    
    # Set up plot symbol.
    shapes <- c(17,16,2,1)
    names(shapes) <- c("distal", "local", "distal_multi", "local_multi")
    x <- dplyr::mutate(x, qtl_type = names(shapes)[1 + local + 2 * (qtl_ct > 1)])
  }

  # Colors
  cols <- c(RColorBrewer::brewer.pal(4, "Dark2"), "#CCCCCC")
  names(cols) <- relabel
  
  switch(type,
         pos_pvalue = {
           p <- ggplot2::ggplot(dplyr::filter(x, x$pvalue <= maxPvalue)) +
             ggplot2::aes(x=pos, y=-log10(pvalue), col=biotype) +
             ggplot2::aes(symbol=symbol, mediation=mediation) +
             ggplot2::facet_grid(~triad) +
             ggplot2::xlab("Position (Mbp)") +
             ggplot2::ylab("-log10 of p-value")
           if(!is.null(pos_t))
             p <- p +
               ggplot2::geom_vline(xintercept = pos_t, col = "darkgrey")
         },
         pvalue_lod = {
           p <- ggplot2::ggplot(dplyr::filter(x, x$pvalue <= maxPvalue)) +
             ggplot2::aes(y=mediation, x=-log10(pvalue), col=biotype) +
             ggplot2::aes(symbol=symbol, position=pos) +
             ggplot2::facet_grid(~triad) +
             ggplot2::geom_hline(yintercept = lod_t, col = "darkgrey") +
             ggplot2::xlab("-log10 of p-value") +
             ggplot2::ylab("Mediation LOD")
         },
         pos_lod = {
           p <- ggplot2::ggplot(x) + 
             ggplot2::aes(y=mediation, x=pos, col=triad) +
             ggplot2::aes(symbol=symbol, pvalue=pvalue, biotype=biotype) +
             ggplot2::geom_hline(yintercept = lod_t, col = "darkgrey") +
             ggplot2::xlab("Position (Mbp)") +
             ggplot2::ylab("Mediation LOD") +
             ggplot2::scale_color_manual(values = cols)
           if(!is.null(pos_t))
             p <- p +
               ggplot2::geom_vline(xintercept = pos_t, col = "darkgrey")
         })
  if(exists("shapes")) {
    p <- p + ggplot2::geom_point(aes(shape = qtl_type), size = 2) +
      ggplot2::aes(chr = chr, qtl_pos = qtl_pos, multi = multi) +
      ggplot2::scale_shape_manual(values = shapes)
  } else {
    p <- p + ggplot2::geom_point(size = 2)
  }
  
  p + ggplot2::ggtitle(main)
}
