sc <- function(bp_cycle, 
               window = TRUE,  # Changed default to TRUE
               pow_exp = 2,
               normalize = FALSE,
               return_spectrum = FALSE) {
  
  # Input validation
  if (length(bp_cycle) < 10) {
    warning("Very short waveform - results may be unreliable")
  }
  
  # 1) Resample to fixed_len points
  fixed_len <- 1024  # Power of 2 for efficient FFT
  x <- signal::resample(bp_cycle, p = fixed_len, q = length(bp_cycle))
  
  # Remove mean (DC component)
  x <- x - mean(x)
  
  # Optional: Remove linear trend as well
  # t <- 1:fixed_len
  # x <- residuals(lm(x ~ t))
  
  if (window) {
    w <- 0.5 - 0.5 * cos(2*pi*(0:(fixed_len-1))/(fixed_len-1))
    x <- x * w
  }
  
  # 2) FFT (fs = fixed_len → frequency axis = harmonic index)
  X      <- fft(x)
  half_N <- fixed_len / 2
  mag    <- Mod(X)[1:half_N]
  freqs  <- 0:(half_N-1)  # harmonic number 0,1,2,...
  
  # Calculate weighted centroid
  weights <- mag^pow_exp
  centroid_harm <- sum(freqs * weights) / sum(weights)
  
  # Centroid in Hz
  if (normalize) {
    centroid_norm <- centroid_harm / (half_N - 1)  # Normalize to [0, 1]
    if (return_spectrum) {
      return(list(
        centroid = centroid_harm,
        centroid_normalized = centroid_norm,
        magnitude = mag,
        frequency = freqs
      ))
    }
    return(centroid_norm)
  }
  
  # Optional: return full spectrum for  diagnostics
  if (return_spectrum) {
    return(list(
      centroid = centroid_harm,
      magnitude = mag,
      frequency = freqs
    ))
  }
  
  return(centroid_harm)
}
