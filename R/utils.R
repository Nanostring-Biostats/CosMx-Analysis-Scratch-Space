#' Calculate row standard deviations without converting to dense matrix.
#' @param x A sparse dgCMatrix
#' @examples 
#'  
#' set.seed(30)
#' sm <- Matrix::sparseMatrix(i=sample(1:30,450,replace=TRUE),
#'                            j=sample(1:20,450,replace=TRUE),x=runif(450))
#' all.equal(rowSDs(sm), apply(sm, 1, sd))
#' 
#' @export
rowSDs <- function(x){
  
  ## std = \sqrt((\sum x_i^2 - n \bar{x}^2)/(n-1))
  if(inherits(x, "sparseMatrix")){
    x2 <- x
    x2@x <- x2@x^2
    x.means <- Matrix::rowMeans(x)
    x2.sums <- Matrix::rowSums(x2)
    nn <- ncol(x)
    rowsds <- sqrt((x2.sums - nn*x.means^2 )/ (nn - 1))
  } else {
    rowsds <- apply(x,1,sd)
  }
  return(rowsds) 
}

message_parallel <- function (...){
  system(sprintf("echo \"%s\"", paste0(..., collapse = "")))
}
