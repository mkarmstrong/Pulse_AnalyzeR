beta_dist <- function(x, alpha, beta = 5, plot = FALSE) {
  
  # start of distribution
  s_d <- which.max(x)
  
  # Generate tau values
  l_o <- length(x) - s_d
  tau <- seq(0, 1, length.out = l_o)
  
  # Calculate PDF values using custom formula
  y1 <- tau^(alpha - 1) * (1 - tau)^(beta - 1)
  
  y2 <- c(rep(0, s_d), y1)
  
  y3 <- (0+(y2-min(y2))*(1-0)/(max(y2)-min(y2)))
  
  # Plot the custom beta distribution curve
  if(isTRUE(plot)) {
    # Plot the gamma distribution curve
    plot(y3, type = "l", col = "blue", lwd = 2,
         main = "Beta Distribution Curve",
         xlab = "X",
         ylab = "Probability Density")
  }
  
  return(y3)
  
}
