root_spline <- function (x, y, y0 = 0, verbose = FALSE) {

  if (is.unsorted(x)) {
    ind <- order(x)
    x <- x[ind]; y <- y[ind]
  }
  
  z <- y - y0
  ## which piecewise linear segment crosses zero?
  k <- which(z[-1] * z[-length(z)] <= 0)
  ## analytical root finding
  xr <- x[k] - z[k] * (x[k + 1] - x[k]) / (z[k + 1] - z[k])
  
  ## make a plot?
  if (verbose) {
    plot(x, y, "l"); abline(h = y0, lty = 2)
    points(xr, rep.int(y0, length(xr)))
  }
  
  ## return roots
  return(xr)
  
}
