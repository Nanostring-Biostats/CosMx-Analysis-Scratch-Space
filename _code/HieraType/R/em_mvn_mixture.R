

#' Fit multivariate gaussian mixture to celltype metagene scores
#' 
#' @param scores a cells x components matrix of metagene scores
#' @param initmu optional components x components matrix of initial mean values
#' @param initsigma optional components length list of components x components covariance matrices used for starting values.
#' @param maxiter maximum iterations before quitting em algorithm.
#' @param rel.eps.convergence if (abs(ll[iter] - ll[iter-1]) / ll[iter]) < rel.eps.convergence, then algorithm is converged.
#' @param prior_prob_level optional vector of posterior probability weights.
#' @param min_pos_score_mean minimum value for 
#' @param max_neg_score_mean maximum value for the mean of off-component metagenes.
#' @param component.variance FALSE is an experimental option.  This should typically be the default TRUE.
#' @param verbose  integer in range of 0-4 controlling verbosity of printed messages.
#'
#'
em_mvn_mixture <- function(scores
                           ,initmu
                           ,initsigma
                           ,maxiter = 1000
                           ,rel.eps.convergence = 1e-6
                           ,component.variance = TRUE
                           ,prior_prob_level = NULL
                           ,verbose = 2
                           ,min_pos_score_mean = 0.1
                           ,max_neg_score_mean = 0
                           ,intercepts = NULL
                           ){

 
  if(!is.null(intercepts)){
    stopifnot(length(intercepts) == ncol(scores))
    for(j in 1:length(intercepts)){
      if(nrow(intercepts[[j]]) != ncol(scores)){
        msg <- paste0("'intercepts' argument should be a list of 1 column matrices,"
                            ,"indicating the minimum score for the 'positive' metagene and maximum score for the other metagenes.")
        msg2 <- paste0("# of rows in intercepts component ", j, " = ", nrow(intercepts[[j]]), ", but should match ncol(scores)=", ncol(scores))
        stop(paste0(msg, "\n", msg2))
      }
      if(!is.null(rownames(intercepts[[j]]))){
        stopifnot(all(rownames(intercepts[[j]]) %in% colnames(scores)))
        intercepts[[j]] <- intercepts[[j]][colnames(scores),]
      }
    }
  } 
  #mu <- do.call(rbind, initmu)
  mu <- initmu 
   
  p <- rep(1/ncol(scores), ncol(scores))
  llvec <- rep(0,maxiter)
  if(!is.null(prior_prob_level)){
    min_prior_prob_level <- min(prior_prob_level[prior_prob_level > 0])
    prior_prob_level <- pmax(prior_prob_level, min_prior_prob_level)
  }
  
  for(iter in 1:maxiter){
    
    if(iter==1){
      ll <- lapply(1:ncol(scores), function(k){
        mvtnorm::dmvnorm(scores
                         #,constraint_mat[k,] * exp(mi[k,])
                         ,mu[k,]
                         ,sigma = initsigma[[k]]
                         ,log = TRUE
        ) + log(p[k])
      })
      ### e-step 
      llmat <- do.call(cbind, ll)
      ppmatdenom <- matrixStats::rowLogSumExps(llmat)
      ppmat <- exp(llmat - ppmatdenom)
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
    for(k in 1:ncol(scores)){
      
      muk <- mstep_muk(muk = mu[k,], scores
                       ,sigk = initsigma[[k]], ppk = ppmat[,k]
                       ,k=k, pos_intercept = min_pos_score_mean #0.02
                       ,neg_intercept = max_neg_score_mean
                       ,intercepts = intercepts[[k]]
                       ) 
      mu[k,] <- muk 
      scoresc <- scale(scores, center = muk, scale=FALSE)
      wt <- sqrt(as.numeric(ppmat[,k]))
      sigk <- crossprod(scoresc * wt) / sum(ppmat[,k])
      initsigma[[k]] <- as.matrix(sigk)
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
    if(verbose > 2){
      message(paste0("iter=",iter,", sigma component:"))
      print(initsigma) 
    } 
    ll <- lapply(1:ncol(scores), function(k){
      mvtnorm::dmvnorm(scores
                       #,constraint_mat[k,] * exp(mi[k,])
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
      llvec[iter] <- sum(ppmatdenom + log(prior_prob_level))
    } else {
      llvec[iter] <- sum(ppmatdenom)
    }
    
    message(paste0("iteration ",iter, ", ll: ",llvec[iter])) 
    
    if(is.na(llvec[iter])) stop("NA likelihood; need to debug")
    if(iter > 1){
      if(llvec[iter] < llvec[iter - 1] && iter > 100) stop("likelihood got worse (?); need to debug")
      rel.eps <- (llvec[iter] - llvec[iter - 1])/llvec[iter]
      if(verbose > 0){
        message(paste0("rel.eps = ", rel.eps))
      }
      if(iter > 10 && abs(rel.eps) < rel.eps.convergence) break
    }
  } 
 
  post_probs <- data.table::data.table(as.matrix(ppmat))
  colnames(post_probs) <- colnames(scores)
  post_probs[,best_score:=do.call(pmax,.SD)]
  post_probs[,best_class:=colnames(post_probs)[which.max(.SD)],by=.I,.SDcols=(1:(ncol(post_probs)-1))]
  post_probs[,cell_ID:=rownames(llmat)] 
  data.table::setcolorder(post_probs, c("cell_ID", "best_class", "best_score"))
 
  rownames(mu) <- paste0(colnames(mu), ".cluster")
  names(p) <- colnames(mu) 
 
  return(list(
    post_probs = post_probs
    ,pihat = p
    ,llmat = llmat
    ,mu_hat = mu                       ## means for each mixture component
    ,sigma_hat = initsigma             ## covariance matrices for each mixture component
    ,loglik_iters = llvec[1:iter]
  ))
}





