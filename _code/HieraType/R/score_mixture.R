
### just the e-step; get posterior probabilities from a fitted model
score_mixture <- function(model, scores, prior_prob_level = NULL){
 
  if(!is.null(prior_prob_level)){
    min_prior_prob_level <- min(prior_prob_level[prior_prob_level > 0])
    prior_prob_level <- pmax(prior_prob_level, min_prior_prob_level)
  }
  
  p <- model$pihat
  ll <- lapply(1:ncol(scores), function(k){
    mvtnorm::dmvnorm(scores
                     ,model$mu_hat[k,]
                     ,sigma = model$sigma_hat[[k]]
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
 
  post_probs <- data.table::data.table(as.matrix(ppmat))
  colnames(post_probs) <- colnames(scores)
  post_probs[,best_score:=do.call(pmax,.SD)]
  post_probs[,best_class:=colnames(post_probs)[which.max(.SD)],by=.I,.SDcols=(1:(ncol(post_probs)-1))]
  post_probs[,cell_ID:=rownames(llmat)] 
  data.table::setcolorder(post_probs, c("cell_ID", "best_class", "best_score"))
  
  return(list(post_probs = post_probs))
}
