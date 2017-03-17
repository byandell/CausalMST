corHat <- function(S.hat) {
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
    tmp <- is.positive.definite(Sig.hat[[i]])
    if (!tmp) {
      Sig.hat[[i]] <- make.positive.definite(Sig.hat[[i]])
    }
    Cor.hat[[i]] <- cov2cor(Sig.hat[[i]])
  }
  Cor.hat
}
