weighted_dicrotic <- function(pw, plot = FALSE) {

  # Get derivatives
  dp1 <- fsg721(pw)
  dp2 <- fsg721(fsg721(pw))
  dp3 <- fsg721(fsg721(fsg721(pw)))
  
  # End index
  end <- length(pw)
  
  # Isolate notch area with 1st derivatives
  nni <- which.min(dp1)
  
  # FIND DICROTIC DEPRESSION ------------------------------------------------
  
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
  hr_bps <- (60/(length(pw)/200))/60
  ft <- which.max(dp2[1:which.max(dp1)])
  tsys <- round((-0.1*hr_bps+0.45)*200) + ft
  data <- data.frame(a = numeric(31), err = numeric(31))
  alpha <- 0
  
  for(i in seq(1.5, 4.5, 0.1)) {
    beta <- beta_dist(pw, i, 5)
    alpha <- i
    if(which.max(beta) >= tsys) {
      break
    }
  }
  
  beta_dis <- beta_dist(pw, alpha, 5)
  weighted_beta <- dp2 * beta_dis
  dic <- which.max(weighted_beta)
  
  # plot(pw, type="l", lwd=2)
  # par(new=T)
  # plot(weighted_beta, type="l",col=2); abline(v=dic, col=2, lty=2)
  # par(new=T)
  # plot(beta_dis, type="l", col=3)
  
  # FIND DICROTIC PEAK ------------------------------------------------------
  
  end3 <- ((end - dic) * .6) + dic # 60% of diastolic duration
  
  # Dicrotic peak from min of 2nd derivative
  # works better for subtle peaks
  # if there are no
  if(sum(dp2[dic:end3] < 0) < 1) {
    dia <- 9999
  } else {
    dia <- which.min(dp2[dic:end3]) + dic
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
  
  # PLOTS -------------------------------------------------------------------
  
  if(isTRUE(plot)) {
    plot(pw, type = "l", lwd=2, ylab="BP (mmHg)")
    abline(v=c(dic, dia), col="grey", lty=3, lwd=2)
    mtext(c("Ed", "P3"), side = 3, at = c(dic,dia))
    par(new = T)
    plot(weighted_beta, type="l", col=2)
    par(new = T)
    plot(beta_dis, type="l", col=3)
  }
  
  return(data.frame(dicrotic_notch = dic, 
                    dicrotic_peak = dia))
  
}
