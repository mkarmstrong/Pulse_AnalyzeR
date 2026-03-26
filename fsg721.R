fsg721 <- function(x, order = 2, smth = 11) {
  sg <- signal::sgolay(p = order, n = smth, m = 1)
  sig <- signal::filter(sg, x)
  return(sig)
}
