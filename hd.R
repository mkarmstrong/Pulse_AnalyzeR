hd <- function(bp_wave, n_harms = 8) {
  
  # Optional: normalize amplitude (doesn't affect HD calculation) but useful for plots
  bp_wave <- (bp_wave - min(bp_wave)) / (max(bp_wave) - min(bp_wave))
  
  # Calculate discrete Fourier transform (DFT)
  fft_coeffs <- fft(bp_wave)
  
  # Calculate power spectrum |Ak|^2
  power_spectrum <- abs(fft_coeffs)^2
  
  # Index 2 is fundamental frequency (k=1)
  fundamental_energy <- power_spectrum[2]
  
  # Sum harmonics from k=2 to k=n_harms (indices 3 to n_harms+1)
  harmonic_energy <- sum(power_spectrum[3:(n_harms + 1)])
  
  # Calculate HD
  hd_value <- harmonic_energy / fundamental_energy
  
  return(hd_value)
}
