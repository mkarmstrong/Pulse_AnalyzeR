low_pass_filt <- function(y, fq = 0.1, do.plot = FALSE) {
  
  if (any(is.na(y))) stop("y contains NA")
  
  # n = a numeric value giving the order of the filter. 
  # Larger numbers create steeper fall off.
  n = 4
  
  if (any(fq>1)) {
    f <- 1/fq
    p <- fq
  } else {
    p <- 1/fq
    f <- fq
  }
  
  # sort f in case it's passed in backwards
  f <- sort(f)
  
  filt <- signal::butter(n = n,
                         W = f * 2,
                         type = "low",
                         plane = "z")
  
  # remove mean
  yAvg <- mean(y)
  y <- y - yAvg
  
  # pad the data to twice the max period
  pad <- max(p) * 2
  ny <- length(y)
  
  # pad the data
  yPad <- c(y[pad:1], y, y[ny:(ny - pad)])
  
  # run the filter
  yFilt <- signal::filtfilt(filt, yPad)
  
  # unpad the filtered data
  yFilt <- yFilt[(pad + 1):(ny + pad)]
  
  # return with mean added back in
  filt.sig <- yFilt + yAvg
  
  if(isTRUE(do.plot)){
    # plot results
    plot(filt.sig,
         type = "l",
         lwd = 2)
  }
  
  # return filtered signal
  return(filt.sig)
  
}
