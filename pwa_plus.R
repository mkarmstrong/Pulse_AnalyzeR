pwa_plus <- function(pw, fs = 200, filt = FALSE, norm = FALSE, verbose = FALSE) {

  # Low pass waveform
  if (isTRUE(filt)) {
    pw <- low_pass_filt(pw, 0.11)
  }

  # Normalize pw amplitude if needed
  if(isTRUE(norm)) {
  pw <- (0+(pw-min(pw))*(1-0)/(max(pw)-min(pw)))
  }
  
  # Calc derivatives
  d1 <- fsg721(pw)
  d2 <- fsg721(d1)
  d3 <- fsg721(d2)
  d4 <- fsg721(d3)
  
  # Remove leading tail, important if pw is averaged using R wave in the ECG.
  rm_lead <- which.max(d2[1:which.max(d1)]) - 1
  if(rm_lead < 1) {rm_lead <- 1}
  pw <- pw[rm_lead:length(pw)]
  ft_i <- 1

  # re-calc derivatives
  d1 <- d1[(rm_lead):(length(d1))]
  d2 <- d2[(rm_lead):(length(d2))]
  d3 <- d3[(rm_lead):(length(d3))]
  d4 <- d4[(rm_lead):(length(d4))]
  
  # plot(fsg721(pw))
  # par(new=T)
  # plot(d1, col = 2)
 
  # Clac time
  sr <- fs
  time <- (0:(length(pw)-1)) / sr

  # Find dicrotic features
  dn_data <- weighted_dicrotic(pw, plot = F) # p @ es
  dic_i <- dn_data$dicrotic_notch              # dic = e
  dia_i <- dn_data$dicrotic_peak
  
  # Find max up slope in systole
  max_d1_i <- which.max(d1)
  
  # Find max (SBP)
  max_pw_i <- which.max(pw)
  
  # Calc sub endocardial variability ratio
  end <- length(pw)
  sevr <- (sum(pw[(dic_i+1):end])/sr) / (sum(pw[round(ft_i):dic_i])/sr)
  
  # Calc harmonic distortion
  har_dist <- hd(pw)
  
  # Calc spectral centroid
  spec_cent <- sc(pw, pow_exp = 1.8)

  # Find a
  a_i <- which.max(d2[1:which.max(d1)]) # basically ft_i
  
  # Find b
  b_i <- which.min(d2[ft_i:max_pw_i])
  
  # Find c
  pks <- find_peaks(d2)
  temp <- which(pks > b_i & pks < dic_i)
  rel_tr_el <- which.max(d2[pks[temp]])
  c_i <- pks[temp[rel_tr_el]]
  if (length(c_i) == 0) {
    trs <- find_peaks(-d3)
    temp <- which(trs > b_i & trs < dic_i)
    rel_tr_el <- which.min(d3[trs[temp]])
    if (length(rel_tr_el) > 0) {
      c_i <- trs[temp[rel_tr_el]]
    }
  }
  
  # Find d
  trs <- find_peaks(-d2)
  possible_trs <- which(trs > c_i & trs < dic_i)
  if (length(possible_trs) > 0) {
    temp <- trs[possible_trs]
    temp_el <- which.min(d2[temp])
    d_i <- temp[temp_el]
  } else {
    d_i <- c_i
  }
  
  # Calc a1 and a2
  a1 = sum(pw[a_i:dic_i])/sr
  a2 = sum(pw[(dic_i+1):end])/sr
  
  # Slopes
  slope_b_c <- unname(coef(lm(d2[b_i:c_i] ~ time[b_i:c_i]))[2])/d2[a_i]
  slope_b_d <- unname(coef(lm(d2[b_i:d_i] ~ time[b_i:d_i]))[2])/d2[a_i]

  # plot(d2[b_i:d_i] ~ time[b_i:d_i])
  # abline(lm(d2[b_i:d_i] ~ time[b_i:d_i]))
  
  # Find pulse wave inflections --------------------------------------------

  # Find early systolic inflection (local max of 2nd derivative per SphygmoCor)
  loc_hold <- which.min(d2[max_d1_i:max_pw_i]) + (max_d1_i - 1)
  p1_in_i <- which.max(d2[loc_hold:max_pw_i]) + (loc_hold - 1)

  # plot(pw,type="b"); abline(v=p1_in_i, col=3)
  # par(new=T)
  # plot(d2,type="l",col="grey"); abline(h=0)

  # Find p2 from 3rd derivative (p2 = max_pw_i or late systolic inflection)
  # When looking for p2 after max_pw_i, using d3 is more robust than using d4.
  p2_i <- which.min(d3[max_pw_i:(dic_i - 5)]) + max_pw_i - 1

  # plot(pw); abline(v=p2_i, col=3)
  # par(new=T)
  # plot(d3, type="l", col="grey"); abline(h=0)
  # abline(v=c(max_pw_i, dic_i-5))

  # Determine waveform type -------------------------------------------------

  # Murgo type is determined from the systolic inflection (ie p1_in_i.
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
  pp  <- (pw[max_pw_i] - pw[ft_i])
  aix <- (ap / pp) * 100

  # diastolic decay ---------------------------------------------------------

  # p & t in diastole
  pd <- pw[dic_i:end]
  td <- 0:(length(pd)-1) / sr
  td <- time[dic_i:end]

  lmod <- lm(pd ~ td)
  tc <- unname(lmod$coefficients["td"])

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
    points(x = c(time[ft_i],
                 time[max_pw_i],
                 time[p2_i],
                 time[dic_i],
                 time[dia_i]),
           y = c(pw[ft_i],
                 pw[max_pw_i],
                 pw[p2_i],
                 pw[dic_i],
                 pw[dia_i]),
           pch = 1,
           col = 2,
           lwd = 2,
           cex = 1.7)
    clip(time[dic_i], time[end], 1000, -1000)
    abline(lmod, col=4, lwd = 2)
    mtext(c("ft", "s", "p2", "dic", "dia"), 
          at = c(time[a_i], time[max_pw_i], time[p2_i], time[dic_i], time[dia_i]),
          side = 3)

    if(aix > 0) {
      points(x = time[p1_in_i],
             y = pw[p1_in_i],
             pch = 1,
             col = 5,
             lwd = 3,
             cex = 1.7)
    }
    
    plot(time, d1,
         type = 'l',
         lwd = 3,
         ylab = "D1 (mmHg/s)",
         xlab = "Time (s)")
    grid(NULL,NULL, lty = 3, col = "lightgrey")
    points(x = c(time[max_d1_i],
                 time[dia_i]),
           y = c(d1[max_d1_i],
                 d1[dia_i]),
           pch = 1,
           col = 2,
           lwd = 2,
           cex = 1.7)
    abline(h=0, lty = 2)
    mtext(c("ms", "dia"), 
          at = c(time[max_d1_i], time[dia_i]),
          side = 3)
    
    plot(time, d2,
         type = 'l',
         lwd = 3,
         ylab = "D2 (mmHg/m^2)",
         xlab = "Time (s)")
    grid(NULL,NULL, lty = 3, col = "lightgrey")
    abline(lm(d2[b_i:c_i] ~ time[b_i:c_i]), col = 2)
    #abline(lm(d2[b_i:d_i] ~ time[b_i:d_i]), col = 2)
    points(x = c(time[a_i],
                 time[b_i],
                 time[c_i],
                 time[d_i],
                 time[dic_i]),
           y = c(d2[a_i],
                 d2[b_i],
                 d2[c_i],
                 d2[d_i],
                 d2[dic_i]),
           pch = 1,
           col = 2,
           lwd = 2,
           cex = 1.7)
    mtext(c("a", "b", "c", "d", "e"), 
          at = c(time[a_i], time[b_i], time[c_i], time[d_i], time[dic_i]),
          side = 3)
    
  }

  # Variables ---------------------------------------------------------------

  # Params per Charlton at DOI:10.1088/1361-6579/aabe6a
  df <- data.frame(
    # PRESSURE
    # timings
    delt_t = time[dia_i] - time[max_pw_i],
    crest_t = time[max_pw_i] - time[a_i], # crest time
    prop_t = (time[max_pw_i] - time[a_i])/time[end],
    sys_t = time[dic_i] - time[a_i],
    dia_t = time[end] - time[dic_i],
    ratio_t = time[max_pw_i]/time[dic_i],
    prop_delt_t = (time[dia_i] - time[max_pw_i])/time[end],
    p1_dia_t = time[dia_i] - time[p1_in_i],
    t_p2_dia = time[dia_i] - time[p2_i],
    ipr = 60/time[end],
    # amplitudes
    ai = (pw[p2_i] - pw[p1_in_i])/pw[max_pw_i],
    ri = pw[dia_i]/pw[p1_in_i],
    ri_p1 = pw[dia_i]/pw[p1_in_i],
    ri_p2 = pw[dia_i]/pw[p2_i],
    ratio_p2_p1 = pw[p2_i]/pw[p1_in_i],
    # areas
    a1,
    a2,
    ipa = (sum(pw[(dic_i+1):end])/sr) / (sum(pw[a_i:dic_i])/sr),
    # DERIVATIVE I
    # amplitudes
    ms = max(d1)/pw[max_pw_i],
    # DERIVATIVE II
    # amplitudes
    b_a = d2[b_i]/d2[a_i],
    c_a = d2[c_i]/d2[a_i],
    d_a = d2[d_i]/d2[a_i],
    e_a = d2[dic_i]/d2[a_i],
    agi = (d2[b_i] - d2[c_i] - d2[d_i] - d2[dic_i])/d2[a_i],
    agi_int = (d2[b_i] - d2[dic_i])/d2[a_i],
    agi_mod = (d2[b_i] - d2[c_i] - d2[d_i])/d2[a_i],
    amb_amp1 = d2[b_i]/pw[p1_in_i],
    # timings
    tb_c = time[c_i] - time[b_i],
    tb_d = time[d_i] - time[b_i],
    # slopes
    slope_b_c,
    slope_b_d,
    # COMBINED
    ipad = (a2/a1) + d_i/a_i,
    k = d2[max_pw_i]/((pw[max_pw_i]-pw[max_d1_i])/pw[max_pw_i]),
    # OTHER ADDITIONS
    har_dist,
    time_con = tc,
    dpdt_slope = max(d1),
    sevr = sevr, # same as ipa above
    ap_mmHg = ap,
    aix = aix,
    p3_sec = time[dia_i],
    p1_alt_sec = time[p1_in_i],
    spec_cent
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
