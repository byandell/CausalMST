qtlMatrix <- function(cross, chr, pos) {
  cross.type <- class(cross)[1]
  le.mar <- length(chr)
  if(le.mar > 0) {
    qtlo <- qtl::makeqtl(cross, chr, pos, what = "prob")

    # Create design matrix (Haley-Knott regression)
    nr <- nrow(qtlo$prob[[1]])
    ng <- length(qtlo$prob)
    if(cross.type == "f2"){
      tmp <- 
        unlist(lapply(qtlo$prob, function(x) cbind(x[,1] - x[,3], x[,2])))
      hkm <- matrix(tmp, nr, 2 * ng)
    }
    if(cross.type == "bc"){
      tmp <- unlist(lapply(qtlo$prob, function(x) x[,1] - x[,2]))
      hkm <- matrix(tmp, nr, ng)
    }
    cbind(rep(1, nr), hkm)
  }
  else {
    n <- nind(cross)
    matrix(1,n,1)    
  }
}
