# fixes the beat so that it spans foot to foot (opposed to ecg R to R)
#' Foot finding based on local max of the 3rd derivative between index 1 and max in 2nd derivative
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

#-----------------------------------------------------------------------------------------------
# create a beta distribution
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

#-----------------------------------------------------------------------------------------------
# calculates the first derivative with a second order polynomial SG filter (per Alun Huges and Kim Parker)
fsg721 <- function(x, order = 2, smth = 11) {
  sg <- signal::sgolay(p = order, n = smth, m = 1)
  sig <- signal::filter(sg, x)
  return(sig)
}

#-----------------------------------------------------------------------------------------------
# applies second order low pass filter
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

#-----------------------------------------------------------------------------------------------
# finds exact point of crossings using linear segments
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

#-----------------------------------------------------------------------------------------------
# finds the dicrotic notch on the pulse waveform using the beta weighted 2nd derivative
weighted_dicrotic <- function(pw, fs = 200, plot = FALSE) {

  # Get derivatives
  dp1 <- fsg721(pw)
  dp2 <- fsg721(fsg721(pw))
  dp3 <- fsg721(fsg721(fsg721(pw)))
  
  # End index
  end <- length(pw)
  
  # Isolate notch area with 1st derivatives
  #nni <- which.min(dp1)
  
  # FIND DICROTIC DEPRESSION 
  
  # # End index without potential perturbation at end diastole
  # end2 <- end * .9
  # 
  # # Dicrotic notch from local dp2 max
  # dic <- which.max(dp2[nni:end2]) + nni - 1
  # 
  # plot(pw, type="l", lwd=2)
  # par(new=T)
  # plot(dp2, type='o',col="grey")
  # abline(v = dic, h = 0)
  
  # per: https://link.springer.com/article/10.1007/s10877-020-00473-3
  # find dicrotic notch
  HR_bps <- fs / length(pw)
  ft <- which.max(dp2[1:which.max(dp1)])
  tsys <- round((-0.1 * HR_bps + 0.45) * fs) + ft
  
  # find alpha analytically
  tPmax <- which.max(pw)
  T_beat <- length(pw)
  tau_wmax <- (tsys - tPmax) / (T_beat - tPmax)
  tau_wmax <- max(0.01, min(0.99, tau_wmax))
  alpha <- (5 * tau_wmax - 2 * tau_wmax + 1) / (1 - tau_wmax)
  alpha <- max(1.5, min(4.5, alpha))
  
  beta_dis <- beta_dist(pw, alpha, 5, plot = F)
  
  # check with standard weighting first
  weighted_beta <- dp2 * beta_dis
  peaks <- which(diff(sign(diff(weighted_beta))) == -2) + 1
  peaks <- peaks[peaks > which.max(pw)]
  
  # plot(pw); abline(v=peaks, lty=2)
  # par(new = T); plot(weighted_beta, type = "l", col=2)
  
  # if the two peaks after sbp are within 75% of each other, sharpen the beta
  # this works better for peripheral waveforms where p2 can get mixed up
  # for central waves the "if" below will always be true (mostly) but it just amplifies the d2 peak associated with dicrotic notch, so not an issue (not robustly tested). 
  # 0.75 is arbitrary, but seems to work quite well
  if (length(peaks) >= 2) {
    top2 <- weighted_beta[peaks][1:2]
    if (top2[1] / top2[2] > 0.75) {
      weighted_beta <- dp2 * beta_dis^2.5
    }
  }
  
  # maybe there is an argument to always use the ^2.5 weighting. 
  # but the original works very well, and changes may/always break stuff
  # if it ain't broke, don't fix it (for now)
  
  dic <- which.max(weighted_beta)
  
  # plot(pw, type="l", lwd=2)
  # par(new=T)
  # plot(weighted_beta, type="l",col=2); abline(v=dic, col=2, lty=2)
  # par(new=T)
  # plot(beta_dis, type="l", col=3)

  
  # FIND DICROTIC PEAK 
  
  end3 <- ((end - dic) * .6) + dic # 60% of diastolic duration
  
  # Dicrotic peak from min of 2nd derivative
  # works better for subtle peaks
  # if there are no
  if(sum(dp2[dic:end3] < 0) < 1) {
    dia <- 9999
  } else {
    dia <- which.min(dp2[dic:end3]) + dic - 1
  }
  
  # plot(pw, type="l", lwd=2)
  # par(new=T)
  # plot(dp2, type='o',col="grey")
  # abline(v = c(dic, dia), col=c(1,2), h = 0)
  
  # Dicrotic peak from 0 crossing of 1st derivative
  # works better for very definable peaks
  # this will over write the value of dia derived above from the 2nd derivative
  if (pw[dia] > pw[dic] & dia != 9999) {
    
    hold <- root_spline(1:length(dp1[(dic):end3]), dp1[(dic):end3])
    dia_hold <- hold[dp1[(dic):end3][hold] > 0]
    # error trap if no 0 crossing is above zero take the max of dp1
    if (length(dia_hold) == 0) {dia_hold <- which.max(dp1[(dic):end3]) - 1}
    
    # if more than 1 peak is found use the highest peak
    if(length(dia_hold) > 1) {
      which_dia_hold <- which.max(pw[dia_hold + dic])
      dia <- dia_hold[which_dia_hold] + dic
    } else {
      dia <- dia_hold + dic
    }
    
  }
  
  # plot(pw, type="l", lwd=2)
  # par(new=T)
  # plot(dp1, type='o',col="grey")
  # abline(v = c(dic, dia), h = 0)
  
  # PLOTS 
  
  if(isTRUE(plot)) {
    plot(pw, type = "l", lwd=2, ylab="BP (mmHg)")
    abline(v=c(dic, dia), col="grey", lty=3, lwd=2)
    mtext(c("Ed", "P3"), side = 3, at = c(dic,dia))
    par(new = T)
    plot(weighted_beta, type="l", col=2)
    par(new = T)
    plot(beta_dis, type="l", col=3)
  }
  
  if (length(dia) == 0) dia <- NA
  if (length(dic) == 0) dic <- NA
  
  return(data.frame(dicrotic_notch = dic, 
                    dicrotic_peak = dia))
  
}

#-----------------------------------------------------------------------------------------------
# harmonic distortion
hd <- function(bp_wave, n_harms = 8) {
  
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

#-----------------------------------------------------------------------------------------------
# Spectral centroid
sc <- function(bp_cycle, window = TRUE, pow_exp = 2, get_spec = FALSE) {
  
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
  
  # Optional: return full spectrum for  diagnostics
  if (get_spec) {
    return(list(
      centroid = centroid_harm,
      magnitude = mag,
      frequency = freqs
    ))
  }
  
  return(centroid_harm)
}

#-----------------------------------------------------------------------------------------------
# find peaks
find_peaks <- function (x, m = 3){
  shape <- diff(sign(diff(x, na.pad = FALSE)))
  pks <- sapply(which(shape < 0), FUN = function(i){
    z <- i - m + 1
    z <- ifelse(z > 0, z, 1)
    w <- i + m + 1
    w <- ifelse(w < length(x), w, length(x))
    if(all(x[c(z : i, (i + 2) : w)] <= x[i + 1])) return(i + 1) else return(numeric(0))
  })
  pks <- unlist(pks)
  pks
}

#-----------------------------------------------------------------------------------------------
# finds intersection of 2 linear regression slopes
lmIntx <- function(fit1, fit2, rnd=2) {
  b1<- fit1$coefficient[1]  #y-int for fit1
  m1<- fit1$coefficient[2]  #slope for fit1
  b2<- fit2$coefficient[1]  #y-int for fit2
  m2<- fit2$coefficient[2]  #slope for fit2
  if(m1==m2 & b1==b2) {print("Lines are identical")
  } else if(m1==m2 & b1 != b2) {print("Lines are parallel")
  } else {
    x <- (b2-b1)/(m1-m2)      #solved general equation for x
    y <- m1*x + b1            #plug in the result
    data.frame(x=round(x, rnd), y=round(y, rnd))
  }
}

#-----------------------------------------------------------------------------------------------
# calculate angles of a triangle
triangle_angles <- function(x1, y1, x2, y2, x3, y3) {
  # Calculate the lengths of the sides using the distance formula
  sideAB <- sqrt((x2 - x1)^2 + (y2 - y1)^2)
  sideBC <- sqrt((x3 - x2)^2 + (y3 - y2)^2)
  sideCA <- sqrt((x1 - x3)^2 + (y1 - y3)^2)
  
  #print(c(sideAB, sideBC, sideCA))
  
  # Calculate the angles using the Law of Cosines
  angleA <- acos((sideBC^2 + sideCA^2 - sideAB^2) / (2 * sideBC * sideCA))
  angleB <- acos((sideCA^2 + sideAB^2 - sideBC^2) / (2 * sideCA * sideAB))
  angleC <- acos((sideAB^2 + sideBC^2 - sideCA^2) / (2 * sideAB * sideBC))
  
  # Convert angles from radians to degrees
  angleA_deg <- angleA * 180 / pi
  angleB_deg <- angleB * 180 / pi
  angleC_deg <- angleC * 180 / pi
  
  #print(c(angleA_deg, angleB_deg, angleC_deg))
  
  # Create a list with angle values
  angles <- list(angleA_deg = angleA_deg, angleB_deg = angleB_deg, angleC_deg = angleC_deg)
  
  return(angles)
}

#-----------------------------------------------------------------------------------------------
# calculat pulse sharpness
psi <- function(pulse, nsample = 500, fs = 200, plot = F) {
  
  pw <- low_pass_filt(pulse, fq = 0.30)
  
  spline_out <- spline(
    x = 1:length(pw),
    y = pw,
    n = nsample)
  
  pw <- spline_out[['y']]

  # create derivatives
  d1 <- fsg721(pw)
  d2 <- fsg721(fsg721(pw))
  
  # get anchor points
  point1 <- which.max(d1)
  d2_zero_cross <- root_spline(1:length(d2), d2)
  point2 <- round(d2_zero_cross[d2_zero_cross > which.max(pw)][1])
  
  # fit left tangent line
  l_x <- (point1-3):(point1+3)
  l_y <- pw[l_x]
  r1 <- abs(cor(l_y, l_x))
  while(r1 > 0.999) {
    l_x <- seq(min(l_x)-1, max(l_x)+1)
    l_y <- pw[l_x]
    r1 <- abs(cor(l_y, l_x))
  }
  mod1 <- lm(l_y ~ l_x)
  m1 <- unname(coef(mod1)[2])
  b1 <- unname(coef(mod1)[1])
  
  # fit right tangent line
  r_x <- (point2-3):(point2+3)
  r_y <- pw[r_x]
  r2 <- abs(cor(r_y, r_x))
  while(r2 > 0.999) {
    r_x <- seq(min(r_x)-1, max(r_x)+1)
    r_y <- pw[r_x]
    r2 <- abs(cor(r_y, r_x))
  }
  mod2 <- lm(r_y ~ r_x)
  m2 <- unname(coef(mod2)[2])
  b2 <- unname(coef(mod2)[1])
  
  # intersection of tangent lines
  x_int <- (b2 - b1) / (m1 - m2)
  y_int <- m1 * x_int + b1
  
  # VH
  VH <- y_int - max(pw)
  
  # end points where tangent lines cross max(pw)
  ep1 <- (max(pw) - b1) / m1
  ep2 <- (max(pw) - b2) / m2
  
  # EPA
  time <- 1:nsample/fs
  tria <- triangle_angles(
    time[round(ep1)], pw[round(ep1)],
    x_int/fs, y_int,
    time[round(ep2)], pw[round(ep2)])
  EPA <- tria$angleB_deg
  
  if(isTRUE(plot)) {
    plot(time, pw, ylim = c(min(pw), y_int), pch=21, bg="grey",
         xlab = "Time (s)", ylab = "Amplitude")
    segments(time[1], b1 + m1*1, time[nsample], b1 + m1*nsample, col=2, lwd=2)
    segments(time[1], b2 + m2*1, time[nsample], b2 + m2*nsample, col=2, lwd=2)
    abline(h = c(max(pw), y_int), col=3, lty=2)
    points(x_int/fs, y_int, pch=19, col=4)
    points(c(time[round(ep1)], time[round(ep2)]),
           c(pw[round(ep1)], pw[round(ep2)]), pch=19, col=4)
    polygon(c(time[round(ep1)], x_int/fs, time[round(ep2)]),
            c(pw[round(ep1)], y_int, pw[round(ep2)]),
            border=4, lwd=2)
  }
  
  # pulse sharpness index
  psi <- 1000 / (EPA * (1 + VH))
  
  df <- data.frame(
    psi = psi,
    epa = EPA,
    vh = VH,
    r1,
    r2
  )
  
  return(df)
}

#-----------------------------------------------------------------------------------------------
# this function puts all of the above together plus some additional calculations
pwa_plus <- function(pw, fs = 200, ecgGated = F, filt = FALSE, norm = TRUE, verbose = FALSE) {

  # Low pass waveform
  if (isTRUE(filt)) {
    pw <- low_pass_filt(pw, 0.11)
  }

  # Remove leading tail and stick it on the end, important if pw is averaged using R wave in the ECG like the nihem does.
  
  if(isTRUE(ecgGated)) {
    pw <- beat_fix_nihem(pw)
  }
  
  # Normalize pw amplitude if needed
  if(isTRUE(norm)) {
    #pw <- (pw - mean(pw)) / sd(pw)             # zero normalization
    pw <- (pw - min(pw)) / (max(pw) - min(pw)) # min0-max1 normalization
  }
  
  # Calc derivatives
  d1 <- fsg721(pw)
  d2 <- fsg721(d1)
  d3 <- fsg721(d2)
  d4 <- fsg721(d3)
  
  # plot(fsg721(pw))
  # par(new=T)
  # plot(d1, col = 2)
 
  # Clac time
  sr <- fs
  time <- (0:(length(pw)-1)) / sr
  end <- length(pw)
  
  # find dicrotic features
  dn_data <- weighted_dicrotic(pw, plot = F) # p @ es
  e_i <- dn_data$dicrotic_notch              # dic = e
  f_i <- dn_data$dicrotic_peak
  
  # find max up slope in systole
  max_d1_i <- which.max(d1)
  
  # find max (SBP)
  max_pw_i <- which.max(pw)
  
  # Calc harmonic distortion
  har_dist <- hd(pw)
  
  # Calc spectral centroid
  spec_cent <- sc(pw, pow_exp = 1.8)

  # calc pulse sharpness
  pulse_sharp <- psi(pw, fs = fs, plot = F)
  p_sharp <- pulse_sharp$psi
  
  # Find a
  a_i <- which.max(d2[1:which.max(d1)]) # p foot
  
  # Find b
  nds <- find_peaks(-d2)
  b_i <- nds[nds > a_i][1]
  #b_i <- which.min(d2[a_i:max_pw_i]) + a_i - 1
  
  # Find c
  # c and e can be confused so just enforce a min distance so c doesn't coincide with e
  # this is because i calculate e with a previous algo and not the find_peaks function
  pks <- find_peaks(d2)
  c_candidates <- pks[pks > b_i & pks < e_i]
  min_ce_separation <- 3
  c_candidates <- c_candidates[c_candidates < (e_i - min_ce_separation)]
  c_i <- c_candidates[which.max(d2[c_candidates])]
  if (length(c_i) == 0) {
    trs <- find_peaks(-d3)
    temp <- which(trs > b_i & trs < e_i)
    rel_tr_el <- which.min(d3[trs[temp]])
    if (length(rel_tr_el) > 0) {
      c_i <- trs[temp[rel_tr_el]]
    }
  }
  
  # Find d
  trs <- find_peaks(-d2)
  possible_trs <- which(trs > c_i & trs < e_i)
  if (length(possible_trs) > 0) {
    temp <- trs[possible_trs]
    temp_el <- which.min(d2[temp])
    d_i <- temp[temp_el]
  } else {
    d_i <- c_i
  }
  
  # Calc a1 and a2
  a1 = sum(pw[a_i:e_i])/sr
  a2 = sum(pw[(e_i+1):end])/sr
  
  # Calc sub endocardial variability ratio (useless when using norm = T)
  sevr <- (a2 / a1)
  
  # Slopes
  slope_b_c <- (d2[b_i] - d2[c_i]) / (time[b_i] - time[c_i])/d2[a_i]
  slope_b_d <- (d2[d_i] - d2[b_i]) / (time[d_i] - time[b_i])/d2[a_i]

  # plot(d2[b_i:d_i] ~ time[b_i:d_i])
  # lines(x = c(time[b_i], time[c_i]), y = c(d2[b_i], d2[c_i]), col = 4, lwd = 2)

  # Find pulse wave inflections --------------------------------------------

  # Find early systolic inflection (local max of 2nd derivative per SphygmoCor)
  loc_hold <- which.min(d2[max_d1_i:max_pw_i]) + (max_d1_i - 1)
  p1_in_i <- which.max(d2[loc_hold:max_pw_i]) + (loc_hold - 1)

  # plot(pw,type="b"); abline(v=p1_in_i, col=3)
  # par(new=T)
  # plot(d2,type="l",col="grey"); abline(h=0)

  # Find p2 from 3rd derivative (p2 = max_pw_i or late systolic inflection)
  # When looking for p2 after max_pw_i, using d3 is more robust than using d4.
  p2_i <- which.min(d3[max_pw_i:(e_i - 5)]) + max_pw_i - 1

  # plot(pw); abline(v=p2_i, col=3)
  # par(new=T)
  # plot(d3, type="l", col="grey"); abline(h=0)
  # abline(v=c(max_pw_i, e_i-5))

  # Determine waveform type -------------------------------------------------

  # If d2 is >= 0, then the early systolic p1 inflection was present.
  if(p1_in_i != max_pw_i) {
    # inflection < max_pw_i
    p2_i <- max_pw_i
  } else {
    # inflection > max_pw_i
    p1_in_i <- max_pw_i
    p1i <- max_pw_i
  }

  # plot(pw); abline(v=c(p1_in_i, p2_i),  col=2:3)
  # par(new=T)
  # plot(d2, type="o", col="grey"); abline(h=0)

  # Calculate augmentation index
  ap  <- (pw[p2_i] - pw[p1_in_i])
  pp  <- (pw[max_pw_i] - pw[a_i])
  aix <- (ap / pp) * 100

  # diastolic decay ---------------------------------------------------------

  # p & t in diastole
  pd_og <- pw[f_i:end]
  td_og <- time[f_i:end]
  
  n <- length(pd_og)
  idx <- round(n * 0.1):round(n * 0.9)
  
  pd <- pd_og[idx]
  td <- td_og[idx]
  
  lmod <- lm(pd ~ td)
  lin_slope <- unname(lmod$coefficients["td"])
  # plot(td, pd, pch = 19, cex = 0.5)
  # lines(td, predict(lmod), col = "red", lwd = 2)
  
  # Plot --------------------------------------------------------------------

  # Plot results
  if (isTRUE(verbose)) {

    par(mfrow=c(3,1),
        mar = c(3.5, 3.5, 1.5, .5),
        mgp = c(2, 1, 0))
    
    plot(time, pw,
         type = 'l',
         lwd = 3,
         ylab = "Pressure (mmHg)",
         xlab = "Time (s)")
    grid(NULL,NULL, lty = 3, col = "lightgrey")
    points(x = c(time[a_i],
                 time[max_pw_i],
                 time[p2_i],
                 time[e_i],
                 time[f_i]),
           y = c(pw[a_i],
                 pw[max_pw_i],
                 pw[p2_i],
                 pw[e_i],
                 pw[f_i]),
           pch = 1,
           col = 2,
           lwd = 2,
           cex = 1.7)
    mtext(c("ft", "s", "p2", "dic", "dia"), 
          at = c(time[a_i], time[max_pw_i], time[p2_i], time[e_i], time[f_i]),
          side = 3)

    if(aix > 0) {
      points(x = time[p1_in_i],
             y = pw[p1_in_i],
             pch = 1,
             col = 5,
             lwd = 3,
             cex = 1.7)
    }
    clip(time[idx[1]] + time[f_i], time[idx[length(idx)]] + time[f_i], 1000, -1000)
    abline(lmod, col=4, lwd = 2)
    
    plot(time, d1,
         type = 'l',
         lwd = 3,
         ylab = "D1 (mmHg/s)",
         xlab = "Time (s)")
    grid(NULL,NULL, lty = 3, col = "lightgrey")
    points(x = c(time[max_d1_i],
                 time[f_i]),
           y = c(d1[max_d1_i],
                 d1[f_i]),
           pch = 1,
           col = 2,
           lwd = 2,
           cex = 1.7)
    abline(h=0, lty = 2)
    mtext(c("ms", "dia"), 
          at = c(time[max_d1_i], time[f_i]),
          side = 3)
    
    plot(time, d2,
         type = 'l',
         lwd = 3,
         ylab = "D2 (mmHg/m^2)",
         xlab = "Time (s)")
    grid(NULL,NULL, lty = 3, col = "lightgrey")
    lines(x = c(time[b_i], time[c_i]), y = c(d2[b_i], d2[c_i]), col = 4, lwd = 2)
    lines(x = c(time[b_i], time[d_i]), y = c(d2[b_i], d2[d_i]), col = 4, lwd = 2)
    points(x = c(time[a_i],
                 time[b_i],
                 time[c_i],
                 time[d_i],
                 time[e_i],
                 time[f_i]),
           y = c(d2[a_i],
                 d2[b_i],
                 d2[c_i],
                 d2[d_i],
                 d2[e_i],
                 d2[f_i]),
           pch = 1,
           col = 2,
           lwd = 2,
           cex = 1.7)
    mtext(c("a", "b", "c", "d", "e", "f"), 
          at = c(time[a_i], time[b_i], time[c_i], time[d_i], time[e_i], time[f_i]),
          side = 3)
    
    par(mfrow=c(1,1))
    
  }

  # Variables ---------------------------------------------------------------

  # Save pulse wave features into a data frame
  df <- data.frame(
    # PRESSURE
    # timings
    delt_t = time[f_i] - time[max_pw_i],
    crest_t = time[max_pw_i] - time[a_i], # crest time
    prop_t = (time[max_pw_i] - time[a_i])/time[end],
    sys_t = time[e_i] - time[a_i],
    dia_t = time[end] - time[e_i],
    ratio_t = time[max_pw_i]/time[e_i],
    prop_delt_t = (time[f_i] - time[max_pw_i])/time[end],
    p1_dia_t = time[f_i] - time[p1_in_i],
    p2_dia_t = time[f_i] - time[p2_i],
    ipr = 60/time[end],
    # amplitudes
    ai = (pw[p2_i] - pw[p1_in_i])/pw[max_pw_i],
    ri = pw[f_i]/pw[p1_in_i],
    ri_p1 = pw[f_i]/pw[p1_in_i],
    ri_p2 = pw[f_i]/pw[p2_i],
    ratio_p2_p1 = pw[p2_i]/pw[p1_in_i],
    ratio_b_p1 = d2[b_i]/pw[p1_in_i],
    # areas
    a1,
    a2,
    ipa = sevr,
    # DERIVATIVE I
    # amplitudes
    ms = max(d1)/pw[max_pw_i],
    # DERIVATIVE II
    # amplitudes
    ba = d2[b_i]/d2[a_i],
    ca = d2[c_i]/d2[a_i],
    da = d2[d_i]/d2[a_i],
    ea = d2[e_i]/d2[a_i],
    agi = (d2[b_i] - d2[c_i] - d2[d_i] - d2[e_i])/d2[a_i],
    agi_int = (d2[b_i] - d2[e_i])/d2[a_i],
    agi_mod = (d2[b_i] - d2[c_i] - d2[d_i])/d2[a_i],
    b_p1 = d2[b_i]/pw[p1_in_i],
    # timings
    bc_t = time[c_i] - time[b_i],
    bd_t = time[d_i] - time[b_i],
    # slopes
    bc_slope = slope_b_c,
    bd_slope = slope_b_d,
    # COMBINED
    ipad = (a2/a1) + d_i/a_i,
    k = d2[max_pw_i]/((pw[max_pw_i]-pw[max_d1_i])/pw[max_pw_i]),
    # OTHER ADDITIONS
    har_dist,
    spec_cent,
    p_sharp,
    dbp_slope = lin_slope,
    dpdt_slope = max(d1),
    p3_t = time[f_i],
    # REDUNDANT INCLUDED FOR COMPLETNESS
    p1_alt_sec = time[p1_in_i], # p1 shoulder (opposed to inflection)
    ap, # augmentation pressure
    aix = aix # augmentation index as a percentage (see ai)
  )

  # round values in df
  num_cols <- unlist(lapply(df, is.numeric)) # Identify numeric cols
  df[num_cols] <-  round(df[num_cols], 3)    # round numeric cols

  # print values to console
  if (isTRUE(verbose)) {
    for (i in 1:length(df)) {
      print(paste0(names(df[i]), ": ", df[1, i]), quote = F)
    }
  }

  return(df)

}

