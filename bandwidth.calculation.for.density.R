bandwidth.calculation.for.density <- function (x) {
    if (length(x) < 2L) 
        stop("need at least 2 data points")
    hi <- sd(x)
    miou <- diff(quantile(as.numeric(x), c(0.25, 0.75), na.rm = FALSE, names = FALSE, 
    type = 6))
    if (!(lo <- min(hi, miou/1.34))) 
        (lo <- hi) || (lo <- abs(x[1L])) || (lo <- 1)
    0.9 * lo * length(x)^(-0.2)
}
