

#' Normalize cell expression in a raw counts matrix by their totalcounts
#' 
#' @param counts_matrix a cells x genes matrix of raw counts to be normalized
#' @param tc optional vector of totalcounts.  Useful if providing a counts_matrix based on a subset of genes.  
#' If not provided, totalcounts per cells is taken to be the rowSums of the counts matrix.
#' 
#' @export
totalcount_norm <- function(counts_matrix, tc = NULL){
  if(is.null(tc)) tc <- Matrix::rowSums(counts_matrix)
  scale.factor <- mean(tc)
  tc[tc==0] <- 1
  return(Matrix::Diagonal(x = scale.factor/tc, names = TRUE) %*% counts_matrix)
}
