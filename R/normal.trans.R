normal.trans <- function (x) 
{
  x <- rank(x, na.last = "keep")
  qnorm(x/(1 + sum(!is.na(x))))
}

