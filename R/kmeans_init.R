

#' function to get starting values for cluster means and covariances
#' @param scores cells x metagenes matrix of scores
#' @param k number of kmeans components
kmeans_init <- function(scores, k=4){
  
  ### initialize celltype cluster means and covariances
  kmod <- stats::kmeans(scores, centers = k)

  initmu <-  
    lapply(1:k, function(kk){
      apply(scores[kmod$cluster==kk,,drop=FALSE], 2, mean)
    }) 
  initmu <- do.call(rbind, initmu)
  
  initsigma <- 
    lapply(1:k, function(kk){
      cov(scores[kmod$cluster==kk,,drop=FALSE])
    })
  
  return(list(initmu = initmu
              ,initsigma = initsigma
              ,kmod = kmod))
}
