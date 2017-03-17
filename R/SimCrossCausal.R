#' Simulate cross object for CMST example
#' 
#' This simulates a \code{cross} object with 3 phenotypes where y1 has a causal effect
#' on both y2 and y3. See vignette for details.
#' + \epsilon$ and $y_3 = 0.5 \, y_1 + \epsilon$). The length of beta controls the number of phenotypes to be simulated.
#' 
#' @param n.ind number of individuals
#' @param len length(s) of chromosome(s)
#' @param n.mar markers per chromosome
#' @param beta causal effect(s) (slope(s)) of y1 on the other phenotypes
#' @param add.eff,dom.eff additive and dominance effects  
#' @param sig.1,sig.2 residual variances for y1 and other responses
#' @param eq.spacing use equal marker spacing if \code{TRUE}
#' @param cross.type cross type (only "bc" and "f2" implemented)
#' @param normalize use normal scores if \code{TRUE}
#' 
#' @export
#' 
SimCrossCausal <- function(n.ind, len, n.mar, beta, add.eff, dom.eff, 
                           sig2.1 = 1, sig2.2 = 1, eq.spacing = FALSE, 
                           cross.type = c("bc", "f2"), normalize = FALSE) {
  n.traits <- length(beta)
  beta <- matrix(rep(beta, each = n.ind), n.ind, n.traits)
  Map <- sim.map(len, n.mar, eq.spacing = eq.spacing, include.x = FALSE)
  Cross <- sim.cross(map = Map, n.ind = n.ind, type = cross.type)
  mygeno <- pull.geno(Cross)
  q <- mygeno[, "D1M51"]
  
  cross.type <- match.arg(cross.type)
  if (cross.type == "bc") {
    add.q <- q - 1.5
    y1 <- add.q * add.eff + rnorm(n.ind, 0, sqrt(sig2.1))
  }
  if (cross.type == "f2") {
    add.q <- q - 2
    dom.q <- (1 + add.q) * (1 - add.q) - 0.5
    y1 <- add.q * add.eff + dom.q * dom.eff + rnorm(n.ind, 0, sqrt(sig2.1))
  }
  y <- beta * y1 + matrix(rnorm(n.ind * n.traits, 0, sqrt(sig2.2)), n.ind, n.traits)
  y <- data.frame(y1, y)
  names(y) <- paste("y", 1 : (n.traits + 1), sep = "")
  if (normalize) {
    apply(y, 2, normal.trans)
  }
  Cross$pheno <- y
  Cross
}
