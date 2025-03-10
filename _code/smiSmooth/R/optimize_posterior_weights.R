
loglikopt_tau_alpha <- function(data, par){
  yi <- data$yi
  xbhat <- data$xbhat
  wts <- data$wts
  eps <- yi - xbhat
  logtau2 <- par[1]
  logalpha <- par[2] 
  tau2 <- exp(logtau2)
  alpha <- exp(logalpha)
  sigma2i <- abs(xbhat) * alpha 
   
  ll <-  
  -0.5 * log(tau2 + sigma2i) -
    0.5 * ((tau2*yi^2 + xbhat^2*sigma2i)/(tau2*sigma2i) - 
             (xbhat*sigma2i + tau2*yi)^2 / (tau2*sigma2i*(tau2 + sigma2i)))
  
  nll <- -1 * sum(ll * wts)
  return(nll)
}


loglikopt_tau_alphav2 <- function(data, par){
  yi <- data$yi
  xbhat <- data$xbhat
  wts <- data$wts
  eps <- yi - xbhat
  logtau2 <- par[1]
  logalpha <- par[2] 
  tau2 <- exp(logtau2)
  alpha <- exp(logalpha)
  tau2 <- abs(xbhat) * tau2
  sigma2i <- alpha
  
  ll <-  
  -0.5 * log(tau2 + sigma2i) -
    0.5 * ((tau2*yi^2 + xbhat^2*sigma2i)/(tau2*sigma2i) - 
             (xbhat*sigma2i + tau2*yi)^2 / (tau2*sigma2i*(tau2 + sigma2i)))
  
  nll <- -1 * sum(ll * wts)
  return(nll)
}



loglikopt_tau_alphav3 <- function(data, par){
  yi <- data$yi
  xbhat <- data$xbhat
  wts <- data$wts
  eps <- yi - xbhat
  logtau2 <- par[1]
  tau2 <- exp(logtau2)
  sigma2i <- tau2 * eps^2 
   
  ll <-  
  -0.5 * log(tau2 + sigma2i) -
    0.5 * ((tau2*yi^2 + xbhat^2*sigma2i)/(tau2*sigma2i) - 
             (xbhat*sigma2i + tau2*yi)^2 / (tau2*sigma2i*(tau2 + sigma2i)))
  
  nll <- -1 * sum(ll * wts)
  return(nll)
}

loglikopt_tau_alphav4 <- function(data, par){
  yi <- data$yi
  xbhat <- data$xbhat
  wts <- data$wts
  eps <- yi - xbhat
  logtau2 <- par[1]
  logsigma2 <- par[2]
  tau2 <- exp(logtau2)
  sigma2i <- exp(logsigma2) * eps^2 
   
  ll <-  
  -0.5 * log(tau2 + sigma2i) -
    0.5 * ((tau2*yi^2 + xbhat^2*sigma2i)/(tau2*sigma2i) - 
             (xbhat*sigma2i + tau2*yi)^2 / (tau2*sigma2i*(tau2 + sigma2i)))
  
  nll <- -1 * sum(ll * wts)
  return(nll)
}

loglikopt_tau_alphav5 <- function(data, par){
  yi <- data$yi
  xbhat <- data$xbhat
  wts <- data$wts
  eps <- yi - xbhat
  logtau2 <- par[1]
  tau2 <- exp(logtau2)
  sigma2i <- mean(eps^2)
   
  ll <-  
  -0.5 * log(tau2 + sigma2i) -
    0.5 * ((tau2*yi^2 + xbhat^2*sigma2i)/(tau2*sigma2i) - 
             (xbhat*sigma2i + tau2*yi)^2 / (tau2*sigma2i*(tau2 + sigma2i)))
  
  nll <- -1 * sum(ll * wts)
  return(nll)
}


add_posteriormean <- function(metagenes, v, wts=NULL, starts = c(log(0.5), log(0.5))){
  
  if(is.null(wts)) wts <- rep(1, length(metagenes[[1]]$y))
  metagenes <- 
  lapply(metagenes, function(xx, version=v, wt = wts){
  
    if(is.null(wt)) wts <- rep(1, length(xx$yhat)) 
    if(version=="v1"){
      wtsfit <- 
      nlminb(start = c(starts), objective = loglikopt_tau_alpha
             ,data = list(yi =xx$y
                          ,xbhat = xx$yhat
                          ,wts = wt)
             )
      # v1 
      tau2hat <- exp(wtsfit$par[1])
      s2hat <- abs(xx$yhat) * exp(wtsfit$par[2])
        
    } 
    if(version=="v2"){
      wtsfit <- 
      nlminb(start = c(starts), objective = loglikopt_tau_alphav2
             ,data = list(yi =xx$y
                          ,xbhat = xx$yhat
                          ,wts = wt)
             )
      # v2
      tau2hat <- exp(wtsfit$par[1]) * abs(xx$yhat)
      s2hat <- exp(wtsfit$par[2])
       
    }
    if(version=="v3"){
      wtsfit <- 
      nlminb(start = c(starts), objective = loglikopt_tau_alphav3
             ,data = list(yi =xx$y
                          ,xbhat = xx$yhat
                          ,wts = wt)
             )
    
      # v3  
      tau2hat <- exp(wtsfit$par[1]) 
      s2hat <- exp(wtsfit$par[1]) * (xx$y - xx$yhat)^2
    }
    
    if(version=="v4"){
      wtsfit <- 
      nlminb(start = c(starts), objective = loglikopt_tau_alphav4
             ,data = list(yi =xx$y
                          ,xbhat = xx$yhat
                          ,wts = wt)
             )
    
      # v4
      tau2hat <- exp(wtsfit$par[1]) 
      s2hat <- exp(wtsfit$par[2]) * (xx$y - xx$yhat)^2
    }
    
    if(version=="v5"){
      wtsfit <- 
      nlminb(start = c(starts), objective = loglikopt_tau_alphav3
             ,data = list(yi =xx$y
                          ,xbhat = xx$yhat
                          ,wts = wt)
             )
    
      # v5 debug
      tau2hat <- exp(wtsfit$par[1]) 
      s2hat <- exp(wtsfit$par[2]) * (xx$y - xx$yhat)^2
    }
    
    if(version=="v6"){
      wtsfit <- 
      nlminb(start = c(starts), objective = loglikopt_tau_alphav5
             ,data = list(yi =xx$y
                          ,xbhat = xx$yhat
                          ,wts = wt)
             )
    
      # v5 debug
      tau2hat <- exp(wtsfit$par[1]) 
      s2hat <- mean((xx$y - xx$yhat)^2)
    }
   
   ypost <- xx$y * tau2hat / (tau2hat + s2hat) + xx$yhat * (s2hat) / (tau2hat + s2hat)
   xx$ypost <- ypost
   xx$ypost_fit <- wtsfit
   return(xx)
    
  })
  return(metagenes) 
}

#