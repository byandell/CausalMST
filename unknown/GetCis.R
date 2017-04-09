##############################################################################
GetCis <- function(x, window = 10) {
  xx <- x[x[, 2] == x[, 4],]
  xx <- xx[abs(xx[, 3] - xx[, 5]) <= window, ]
  index <- match(xx[, 1], x[, 1])
  list(cis.reg = xx, cis.index = index) 
}
