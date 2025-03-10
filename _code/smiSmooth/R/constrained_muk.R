

#' Objective function for M-step for 'mu' where one 
#' component mean has a minimum 'positive intercept'
#' and remaining means have maximum 'negative intercept'
#' 
#'  
constrained_muk <- function(param, scores, cholsiginv, ppk, K, Xk, intercepts){

  K <- ncol(scores)
  beta <- matrix(exp(param), ncol=1,nrow=K)
  
  muk <- Xk %*% beta + intercepts
 
  scoresc <- scale(scores, center = muk, scale = FALSE) 
  objective <- Matrix::Diagonal(x=sqrt(ppk)) %*% scoresc  %*% cholsiginv
  objective <- sum(Matrix::colSums(objective^2))

  return(objective)
}


mstep_muk <- function(muk, scores, sigk, ppk, k, pos_intercept = 0.05, neg_intercept=0, intercepts = NULL){
  
  K <- ncol(scores)
  
  if(is.null(intercepts)){
    intercepts <- matrix(neg_intercept, ncol = 1, nrow = K) 
    intercepts[k,] <- pos_intercept
  }
    
  Xk <- diag(x = rep(-1, K))      ## constraints
  Xk[k,k] <- 1
  
  #muk <- Xk %*% beta + intercepts
  siginv <- solve(sigk)
  cholsiginv <- t(chol((siginv)))
  
  
  mukfit <- 
    nlminb(start = c(log(abs(muk))), objective = constrained_muk
           ,scores=scores, cholsiginv=cholsiginv, ppk = ppk, K=K, Xk = Xk, intercepts = intercepts
    )
  
  muknew <- c(Xk %*% exp(mukfit$par) + intercepts)
  names(muknew) <- colnames(scores)
  
  return(c(muknew))
}

#mukfit2 <- 
#  nlminb(start = c(log(abs(muk))), objective = constrained_muk
#         ,scores=scores, sigk=initsigma[[k]], ppk = prior_prob_level * ppmat[,k], k=k, pos_intercept = 0.05
#  )











