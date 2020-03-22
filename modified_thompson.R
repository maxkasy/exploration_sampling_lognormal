#D vector based on shares
ProportionalAssignment = function(Shares, Nt) {
  k = length(Shares)
  n_floor = floor(Nt * Shares) #number of units assigned to each treatment, rounded down
  remainder = Nt * Shares - n_floor
  units_remaining = Nt - sum(n_floor) #remaining number of units to be assigned
  
  if (units_remaining > 0) {
    Dt = c(rep(1:k, n_floor),
           # assigning each treatment acoording to rounded down average count
           sample(
             1:k,
             size = units_remaining,
             replace = T,
             prob = remainder
           )) #remaining units assigned randomly from remain replicate assignments
  } else {
    Dt = rep(1:k, n_floor)
  }
  sample(Dt)
}


DtchoiceThompson = function(Y, D, #outcomes and treatments thus far
                            k, #number of treatments
                            C, #vector of treatment cost
                            Nt) { # number of observations for period t
  
  # setting prior parameters
  mu_0 = 1
  kappa_0 = 1
  nu_0 = 1
  sigma2_0 = 1
  
  # calculate summary statistics for posterior updating
  log_Y = log(Y)
  nn = tabulate(D)
  bar_log_Y = tapply(log_Y, D, mean, default = 0)
  var_log_Y = tapply(log_Y, D, var, default = 0)
  
  # calculate posterior parameters
  # cf Gelman et al, p67ff
  mu_t = mu_0 + (nn/(kappa_0 + nn)) * (bar_log_Y - mu_0)
  kappa_t = kappa_0 + nn
  nu_t = nu_0 + nn
  sigma2_t = (1/nu_t) * (nu_0 * sigma2_0 +
                        (nn-1)*var_log_Y + 
                        ((kappa_0*nn) /(kappa_t)) * (bar_log_Y - mu_0)^2)
  
  # sample from posterior and determine optimal treatment
  Dt = rep(0, Nt)
  
  for (i  in 1:Nt) {
    sig2 = sigma2_t * nu_t / rchisq(n=k, df=nu_t)  # draw sig2 from scaled inverse chi2 posterior
    mu = rnorm(n=k,mu_t, sig2 / kappa_t)  # draw mu from normal posterior given sig2
    thetadraw = exp(mu + sig2 / 2) # calculate mean of the exponential of N(mu, sig2)
  
    Dt[i] = which.max(thetadraw - C)
  }
  
  factor(Dt, levels = 1:k)
}



DtchoiceThompsonProbabilities = function(Y, D, #outcomes and treatments thus far
                                         k, #number of treatments
                                         C = rep(0, k), #vector of treatment cost
                                         RR = 50000) {
  #number of replication draws
  
  # Repeat Thompson sampling RR times
  DtRR = DtchoiceThompson(Y, D, k, C, RR)
  P_Dt = table(DtRR) / RR #average count for each treatment value and covariate value, replicated sample
  
  P_Dt = as_tibble(matrix(P_Dt, 1, k))
  colnames(P_Dt) = paste(1:k)
  P_Dt
}


DtchoiceThompson_modified = function(alpha) {
  if (max(alpha) < 1) {
    Shares = alpha * (1 - alpha) # Based on stationary distribution for modified Thompson
    Shares = Shares / sum(Shares)
  } else {
    Shares = alpha
  }
  
  return(Shares)
}
