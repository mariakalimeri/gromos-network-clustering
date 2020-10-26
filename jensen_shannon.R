jensen_shannon <- function(p1, p2) {
  # This function assumes that p1 and p2 have the same x
  p <- 0.5 * (p1 + p2)
  js <- 0.5 * sum(p1 * log(p1 / p)) + 0.5 * sum(p2 * log(p2 / p))
  # The following part is necessary in case p1 and p2 are the same. In such
  # a case, it is possible due to rounding errors, that the log(p1/p) or
  # log(02/p) is almost zero (e-15) but from the zero side and that will
  # mess up the sqrt afterwards.
  # Attention: It is the square root of js that is a metric
  if (!is.nan(js) && js < 0 && abs(js) < 1e-10) {
    js <- -js
  }
  return(js)
}
