oldShat <- function(vec.logLik) {
  n.ind <- nrow(vec.logLik)
  vec.LR <- matrix(NA, n.ind, 6, 
                   dimnames = list(NULL, c("12", "13", "14", "23", "24", "34")))
  ip <- matrix(c(1, 1, 1, 2, 2, 3, 2, 3, 4, 3, 4, 4), 6, 2)
  for (i in 1 : 6) {
    vec.LR[, i] <- vec.logLik[, ip[i, 1]] - vec.logLik[, ip[i, 2]]
  }
  S.hat <- (1 - 1 / n.ind) * cov(vec.LR)
  dimnames(S.hat) <- list(dimnames(vec.LR)[[2]], dimnames(vec.LR)[[2]])
  S.hat
}
oldCalcZ <- function(S.hat, ICs, n.ind) {
  Z.ic <- matrix(NA, 4, 4)
  ii <- 1
  for (i in 1 : 3) {
    for (j in (i + 1) : 4) {
      LRt <- -0.5 * (ICs[i] - ICs[j])
      Z.ic[i, j] <- LRt / sqrt(S.hat[ii, ii] * n.ind)
      ii <- ii + 1
    }
  }
  Z.ic
}
oldCorHat <- function(S.hat) {
  Sig.hat <- vector(mode = "list", length = 4) 
  signs.2 <- matrix(c(1, -1, -1, -1, 1, 1, -1, 1, 1), 3, 3)
  signs.3 <- matrix(c(1, 1, -1, 1, 1, -1, -1, -1, 1), 3, 3)
  # 1: 12.13.14
  Sig.hat[[1]] <- S.hat[c(1, 2, 3), c(1, 2, 3)]
  # 2: 21.23.24
  Sig.hat[[2]] <- S.hat[c(1, 4, 5), c(1, 4, 5)] * signs.2
  # 3: 31.32.34
  Sig.hat[[3]] <- S.hat[c(2, 4, 6), c(2, 4, 6)] * signs.3
  # 4: 41.42.43
  Sig.hat[[4]] <- S.hat[c(3, 5, 6), c(3, 5, 6)]
  Cor.hat <- vector(mode = "list", length = 4) 
  for (i in 1 : 4) {
    tmp <- corpcor::is.positive.definite(Sig.hat[[i]])
    if (!tmp) {
      Sig.hat[[i]] <- corpcor::make.positive.definite(Sig.hat[[i]])
    }
    Cor.hat[[i]] <- cov2cor(Sig.hat[[i]])
  }
  Cor.hat
}
oldParCMST <- function(Z) {
  pv <- matrix(NA, 4, 4)
  for (i in 1 : 3) {
    for(j in (i + 1) : 4) {
      pv[i, j] <- pnorm(Z[i, j], lower.tail = FALSE)
    }
  }
  pval.1 <- max(pv[1, 2], pv[1, 3], pv[1, 4])
  pval.2 <- max(1 - pv[1, 2], pv[2, 3], pv[2, 4])
  pval.3 <- max(1 - pv[1, 3], 1 - pv[2, 3], pv[3, 4])
  pval.4 <- max(1 - pv[1, 4], 1 - pv[2, 4], 1 - pv[3, 4])
  c(pval.1, pval.2, pval.3, pval.4)
}
oldNonparCMST <- function(penalty, n, k, vec.logLik) {
  # Computes the non-parametric CMST
  # k: vector of length 4 with the model dimensions
  # vec.logLik: matrix of loglik scores (n by 4)
  kM <- matrix(rep(k, each = n), n, 4)
  if (penalty == "bic") {
    vec.penal <- vec.logLik - 0.5 * kM * log(n)/n
  }
  if (penalty == "aic") {
    vec.penal <- vec.logLik - kM/n
  }
  pv <- matrix(NA, 4, 4)
  for (i in 1 : 3) {
    for (j in (i + 1) : 4) {
      vec.ratio <- vec.penal[, i] - vec.penal[, j]
      counts <- sum(vec.ratio > 0)
      nn <- n - sum(vec.ratio == 0)
      pv[i, j] <- pbinom(counts - 1, nn, 0.5, lower.tail = FALSE)
      pv[j, i] <- pbinom(counts, nn, 0.5, lower.tail = TRUE)
    }
  }
  pval.1 <- max(pv[1, 2], pv[1, 3], pv[1, 4])
  pval.2 <- max(pv[2, 1], pv[2, 3], pv[2, 4])
  pval.3 <- max(pv[3, 1], pv[3, 2], pv[3, 4])
  pval.4 <- max(pv[4, 1], pv[4, 2], pv[4, 3])
  c(pval.1, pval.2, pval.3, pval.4)
}
oldJointCMST <- function(object, Cor.hat) {
  
  z <- min(Z[1, 2], Z[1, 3], Z[1, 4])
  pval.1 <- 1 - mnormt::pmnorm(c(z, z, z), c(0, 0, 0), Cor.hat[[1]])
  z <- min(- Z[1, 2], Z[2, 3], Z[2, 4])
  pval.2 <- 1 - mnormt::pmnorm(c(z, z, z), c(0, 0, 0), Cor.hat[[2]])
  z <- min(-Z[1, 3], -Z[2, 3], Z[3, 4])
  pval.3 <- 1 - mnormt::pmnorm(c(z, z, z), c(0, 0, 0), Cor.hat[[3]])
  z <- min(-Z[1, 4], -Z[2, 4], -Z[3, 4])
  pval.4 <- 1 - mnormt::pmnorm(c(z, z, z), c(0, 0, 0), Cor.hat[[4]])
  c(pval.1, pval.2, pval.3, pval.4)
}

