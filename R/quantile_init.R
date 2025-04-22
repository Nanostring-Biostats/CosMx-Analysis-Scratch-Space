

#' function to get starting values for cluster means and covariances
#' @param scores cells x metagenes matrix of scores
#' @param q upper quantile used to get starting means/covariances for initializing clusters
#' @param minqu set the score to this value if the 'q'th quantile is below it.
quantile_init <- function(scores, q=0.9, minqu = 0.05){
  
  ### initialize celltype cluster means and covariances
  qu <- apply(scores, 2, quantile, q)
  qu[qu < minqu] <- minqu
  
  initmu <- 
    lapply(1:ncol(scores), function(xx){
      apply(scores[scores[,xx] > qu[xx],], 2, mean)
    })
  initmu <- do.call(rbind, initmu)
  
  initsigma <- 
    lapply(1:ncol(scores), function(xx){
      cov(scores[scores[,xx] > qu[xx],])
    })
  
  return(list(initmu = initmu
              ,initsigma = initsigma))
}