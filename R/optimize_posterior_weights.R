
loglikopt_tau_alpha <- function(data, par){
  yi <- data$yi
  xbhat <- data$xbhat
  wts <- data$wts
  tau2    <- exp(par[1])
  alpha   <- exp(par[2])
  sigma2i <- abs(xbhat) * alpha
  # 2026-08: Simplified form of the marginal N(xbhat, tau2 + sigma2i) log-likelihood.
  # Algebraically identical to the original expanded form when sigma2i > 0, but
  # avoids dividing by sigma2i, which is 0 whenever xbhat == 0 (NNLS exact zeros).
  ll  <- -0.5 * log(tau2 + sigma2i) - 0.5 * (yi - xbhat)^2 / (tau2 + sigma2i)
  nll <- -1 * sum(ll * wts)
  return(nll)
}


add_posteriormean <- function(metagenes, v, wts=NULL, starts = c(log(0.5), log(0.5))){
  
  if(is.null(wts)) wts <- rep(1, length(metagenes[[1]][[1]]$y))
  metagenes <- 
  lapply(metagenes, function(mm, version=v, wt = wts){
 
    if(is.null(wt)) wts <- rep(1, length(mm[[1]]$yhat))
    for(ii in 1:length(mm)){
      xx <- mm[[ii]]
      wtsfit <- 
      nlminb(start = c(starts), objective = loglikopt_tau_alpha
             ,data = list(yi =xx$y
                          ,xbhat = xx$yhat
                          ,wts = wt)
             )
      tau2hat <- exp(wtsfit$par[1])
      s2hat <- abs(xx$yhat) * exp(wtsfit$par[2])
      ypost <- xx$y * tau2hat / (tau2hat + s2hat) + xx$yhat * (s2hat) / (tau2hat + s2hat)
      xx$ypost <- ypost
      xx$ypost_fit <- wtsfit
      mm[[ii]] <- xx
    }
    return(mm)
  })
  return(metagenes) 
}

#
