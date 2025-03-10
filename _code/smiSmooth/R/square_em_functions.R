


### params <- c(l1_model_overfit$pihat
##             ,as.vector(t(l1_model_overfit$mu_hat))
##             ,unlist(lapply(l1_model_overfit$sigma_hat, as.vector)))

deconstruct_params <- function(params, K, nc){
  p <- params[1:K]
  nmupar <- K*nc
  mu <- params[(K+1):(K + nmupar)]
  mu <- matrix(mu, nrow = K, ncol = nc, byrow = TRUE)
  sigma_idx_seq <- seq(nmupar + K + 1, length(params), by=nc * nc) 
  initsigma <- vector(mode = 'list',length=K)
  for(jj in 1:(length(sigma_idx_seq))){
    initsigma[[jj]] <- matrix(params[sigma_idx_seq[jj]:(sigma_idx_seq[jj]+nc*nc-1)],nrow=ncol(scores), ncol=ncol(scores)) 
  }
  return(list(initsigma = initsigma, mu = mu, p = p))
}

construct_params <- function(initsigma, mu, p){
  params <- c(p, as.vector(t(mu)), unlist(lapply(initsigma, as.vector)))
  return(params)
}

square_gmm.em <- function(params, scores, K, ppmat = NULL, prior_prob_level = NULL){
  parms <- deconstruct_params(params, K, nc = ncol(scores))
  p <- parms$p
  initsigma <- parms$initsigma
  mu <- parms$mu
 
   
  #ppmat <- get("ppmat", envir=parent.frame()) 
  #if(is.null(ppmat)){
  ll <- lapply(1:K, function(k){
    mvtnorm::dmvnorm(scores
                     ,mu[k,]
                     ,sigma = initsigma[[k]]
                     ,log = TRUE
    ) + log(p[k])
  })
  
  ### e-step 
  llmat <- do.call(cbind, ll)  
  ppmatdenom <- apply(llmat, 1, matrixStats::logSumExp)
  ppmatnum <- exp(llmat)
  ppmat <- Matrix::Diagonal(x=1/exp(ppmatdenom),names=TRUE)%*% ppmatnum  
  isna <- apply(ppmat, 1, function(x) sum(is.na(x) | is.infinite(x)))
  if(sum(isna > 0) > 0){
    browser()
    whichisna <- which(isna > 0)
    maxll <- apply(llmat[whichisna,,drop=FALSE], 1, which.max)
    for(j in 1:length(whichisna)){
      ppmat[whichisna[j],] <- 0
      ppmat[whichisna[j],maxll[j]] <- 1
    }
  } 
  if(!is.null(prior_prob_level)){
    ppmat <- (Matrix::Diagonal(x=prior_prob_level,names=TRUE) %*% ppmat)
  } 
  #}
  
  browser()
  ### prior update
  p <- Matrix::colSums(ppmat)
  p <- p / sum(p)
  
  for(k in 1:K){
    muk <- ppmat[,k] %*% scores / sum(ppmat[,k])
    mu[k,] <- muk 
    scoresc <- scale(scores, center = muk,scale=FALSE)
    sigk <- (Matrix::t(scoresc) %*% Matrix::Diagonal(x=ppmat[,k]) %*% scoresc)
    sigk <- sigk/(sum(ppmat[,k]))
    initsigma[[k]] <- as.matrix(sigk)
  }
  params <- construct_params(initsigma, mu, p)
  return(params) 
}

square_gmm.loglik <- function(params, scores, prior_prob_level, K, ppmat){
  parms <- deconstruct_params(params, K, nc = ncol(scores))
  p <- parms$p
  initsigma <- parms$initsigma
  mu <- parms$mu
  ll <- lapply(1:K, function(k){
    mvtnorm::dmvnorm(scores
                     ,mu[k,]
                     ,sigma = initsigma[[k]]
                     ,log = TRUE
    ) + log(p[k])
  })
  
  ### e-step 
  llmat <- do.call(cbind, ll)  
  ppmatdenom <- apply(llmat, 1, matrixStats::logSumExp)
  if(!is.null(prior_prob_level)){
    llval <- sum(ppmatdenom + log(prior_prob_level))
  } else {
    llval <- sum(ppmatdenom)
  }
  #assign("ppmat", ppmat, envir = parent.frame())
  smiDE:::message_parallel(llval)
  return(-1*llval) 
}




