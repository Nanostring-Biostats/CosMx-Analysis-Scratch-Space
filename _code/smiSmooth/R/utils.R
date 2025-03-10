

#' 
#' @export
totalcount_norm <- function(counts_matrix, tc = NULL){
  if(is.null(tc)) tc <- Matrix::rowSums(counts_matrix)
  scale.factor <- mean(tc)
  tc[tc==0] <- 1
  return(Matrix::Diagonal(x = scale.factor/tc, names = TRUE) %*% counts_matrix)
}
