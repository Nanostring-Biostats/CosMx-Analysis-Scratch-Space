







#' Function for fitting multi-component gaussian mixture models
em_mvn_mixture_overfit <- function(scores
                           ,ncomponents
                           ,initmu = NULL
                           ,initsigma=NULL
                           ,maxiter = 1000
                           ,rel.eps.convergence = 1e-6
                           ,component.variance = TRUE
                           ,prior_prob_level = NULL
                           ,verbose = 0
                           ,min_pos_score_mean = 0.1
                           ,max_neg_score_mean = 0
                           ){


  
  mu <- initmu
  K <- nrow(mu)
  p <- rep(1/K, K)
  llvec <- rep(0,maxiter)
  if(!is.null(prior_prob_level)){
    min_prior_prob_level <- min(prior_prob_level[prior_prob_level > 0])
    prior_prob_level <- pmax(prior_prob_level, min_prior_prob_level)
  }
  
  for(iter in 1:maxiter){
    if(iter==1){
      ll <- lapply(1:K, function(k){
        mvtnorm::dmvnorm(scores
                         ,mu[k,]
                         ,sigma = initsigma[[k]]
                         ,log = TRUE
        ) + log(p[k])
      })
      ### e-step 
      llmat <- do.call(cbind, ll)
      ppmatdenom <- matrixStats::rowLogSumExps(llmat)
      ppmat <- exp(llmat - ppmatdenom)
      isna <- rowSums(is.na(ppmat) | is.infinite(ppmat))
      if(sum(isna > 0) > 0){
        whichisna <- which(isna > 0)
        maxll <- apply(llmat[whichisna,,drop=FALSE], 1, which.max)
        for(j in 1:length(whichisna)){
          ppmat[whichisna[j],] <- 0
          ppmat[whichisna[j],maxll[j]] <- 1
        }
      } 
      if(!is.null(prior_prob_level)){
        ppmat <- ppmat * as.numeric(prior_prob_level)
      }
    }
    
    ### maximization step
    ### prior update
    p <- Matrix::colSums(ppmat)
    p <- p / sum(p)
    
    ### means update 
    ### sigmas update
    for(k in 1:K){
      ppdenom <- sum(ppmat[,k])
      if(ppdenom > 0){
        muk <- ppmat[,k] %*% scores / ppdenom
        mu[k,] <- muk 
        scoresc <- scale(scores, center = muk, scale=FALSE)
        wt <- sqrt(as.numeric(ppmat[,k]))
        sigk <- crossprod(scoresc * wt) / ppdenom
        initsigma[[k]] <- as.matrix(sigk)
      }
    } 
  
    if(!component.variance){
      combined_sigma <- Reduce('+' , initsigma) * (1/ncol(scoresc))
      initsigma <- lapply(1:length(initsigma), function(x){combined_sigma})
    } 
    if(verbose > 1){
      message(paste0("iter=",iter,", mu component:"))
      #print(exp(mi) * constraint_mat) 
      print(mu) 
    } 
    if(verbose > 1){
      message(paste0("iter=",iter,", sigma component:"))
      print(initsigma) 
    } 
    ll <- lapply(1:K, function(k){
      mvtnorm::dmvnorm(scores
                       ,mu[k,]
                       ,sigma = initsigma[[k]]
                       ,log = TRUE
      ) + log(p[k])
    })
    
    ### e-step 
    llmat <- do.call(cbind, ll)
    ppmatdenom <- matrixStats::rowLogSumExps(llmat)
    ppmat <- exp(llmat - ppmatdenom)

    isinfll <- rowSums(is.infinite(llmat) & llmat > 0)
    if(sum(isinfll > 0) > 0){
      whichisinfll <- which(isinfll > 0)
      impdenom <- max(ppmatdenom[-c(whichisinfll)])
      maxll <- apply(llmat[whichisinfll,,drop=FALSE], 1, which.max)
      for(j in 1:length(whichisinfll)){
        ppmatdenom[whichisinfll[j]] <- impdenom
        ppmat[whichisinfll[j],] <- 0
        ppmat[whichisinfll[j],maxll[j]] <- 1
      }
    } 
    
    isna <- rowSums(is.na(ppmat) | is.infinite(ppmat))
    if(sum(isna > 0) > 0){
      whichisna <- which(isna > 0)
      maxll <- apply(llmat[whichisna,,drop=FALSE], 1, which.max)
      for(j in 1:length(whichisna)){
        ppmat[whichisna[j],] <- 0
        ppmat[whichisna[j],maxll[j]] <- 1
      }
    } 
    if(!is.null(prior_prob_level)){
      ppmat <- ppmat * as.numeric(prior_prob_level)
      llvec[iter] <- sum(ppmatdenom + log(prior_prob_level))
    } else {
      llvec[iter] <- sum(ppmatdenom)
    }
   
    if(verbose > 0){
      message(paste0("iteration ",iter, ", ll: ",llvec[iter])) 
    } 
    
    if(is.na(llvec[iter])) stop("NA likelihood; need to debug")
    if(iter > 1){
      rel.eps <- (llvec[iter] - llvec[iter - 1])/llvec[iter]
      if(llvec[iter]==-Inf) break
      if(verbose > 0) message(paste0("rel.eps = ", rel.eps))
      if(iter > 10 && abs(rel.eps) < rel.eps.convergence) break
    }
  } 
  
  post_probs <- data.table::data.table(as.matrix(ppmat))
  rownames(mu) <- paste0("k",1:nrow(mu))
  colnames(post_probs) <- rownames(mu)
  
  post_probs[,best_score:=do.call(pmax,.SD)]
  post_probs[,best_class:=colnames(post_probs)[apply(.SD,1,which.max)],.SDcols=1:(ncol(post_probs)-1)]
#  post_probs[,best_class:=colnames(post_probs)[which.max(.SD)],by=.I,.SDcols=(1:(ncol(post_probs)-1))]
  post_probs[,cell_ID:=rownames(llmat)] 
  data.table::setcolorder(post_probs, c("cell_ID", "best_class", "best_score"))
 
  names(p) <- rownames(mu) 
 
  return(list(
    post_probs = post_probs
    ,pihat = p
    ,llmat = llmat
    ,mu_hat = mu                       ## means for each mixture component
    ,sigma_hat = initsigma             ## covariance matrices for each mixture component
    ,loglik_iters = llvec[1:iter]
  ))
}





