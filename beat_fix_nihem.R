beat_fix_nihem <- function(x, plot = FALSE) {
  
  p <- x
  pw <- low_pass_filt(p, 0.08)
  
  d1 <- fsg721(pw)
  d2 <- fsg721(fsg721(pw))
  d3 <- fsg721(fsg721(fsg721(pw)))
  
  #plot(pw,  type="l"); par(new=T); plot(d2, type="l", col=2); abline(v=which.max(d3[1:which.max(d2)]))
  
  foot_idx <- which.max(d3[1:which.max(d2)])
  
  #plot(pw[1:15]); abline(v=foot_idx)
  
  x1 <- p[1:foot_idx]
  x2 <- p[(foot_idx):length(pw)]
  pw2 <- c(x2, x1)
  pw2 <- low_pass_filt(pw2)
  
  #plot(pw2)

  if(isTRUE(plot)) {
    plot(pw, type="l", lwd=2, col=2)
    lines(p, col=1, lwd=2)
    lines(pw2, col=3, lwd=2)
    abline(v = foot_idx)
  }
  
  return(pw2)
  
}
