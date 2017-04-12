#' @export
#'
wilcIUCMST <- function(models, flavor = c("B","A")) {

  # Penalize individual log likelihood ratios
  flavor <- match.arg(flavor)
  d <- dim(models$indLR)
  pen <- penalty(d[1], flavor) / (2 * d[1])
  indLR <- as.data.frame(t(t(models$indLR) - models$df * pen))

  # Outer differences of models (columns) of indLR
  tmpfn <- function(d) {
    nz <- (d != 0)
    r <- rank(abs(d[nz]))
    s <- sign(d[nz])
    c(W = sum(r * s),
      nz = sum(nz))
  }
  # var = nz * (nz + 1) * (2 * nz + 1) / 6

  LR <- dplyr::bind_cols(
    purrr::map(indLR,
               function(x,y) x - y,
               indLR))

  names(LR) <- paste(rep(names(models$indLR), each = ncol(models$indLR)),
                     names(LR),
                     sep = ":")
  # Reduce to unique model comparisons Gj/Gk with j<k.
  LR <- LR[, rep(seq(d[2]), each = d[2]) < rep(seq(d[2]), d[2])]

  LR <- t(data.frame(purrr::map(LR, tmpfn),
                      check.names = FALSE))

  ## Organize as left_right, use pos to get pvalue
  LR2 <- dplyr::mutate(
    dplyr::rename(
      dplyr::mutate(
        left_right(unlist(LR[,"W"])),
        nz = rep(as.vector(LR[,"nz"]), 2),
        v = nz * (nz + 1) * (2 * nz + 1) / 6),
      W = Z),
    pv = pnorm(W, 0, sqrt(v), lower.tail = FALSE))

  # Compare reference model with all others and get max pvalue.
  dplyr::mutate(
    dplyr::ungroup(
      dplyr::summarize(
        dplyr::group_by(
          dplyr::mutate(LR2,
                        ref = factor(ref, unique(ref))),
          ref),
        alt = alt[which.max(pv)][1],
        pv = max(pv))),
    ref = as.character(ref))
}