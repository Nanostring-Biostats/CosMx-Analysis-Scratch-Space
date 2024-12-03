#' Convert dense and dgT matrices to dgC matrices
#'
#' @param x Matrix, dgTMatrix format or dense format
#'
#' @import Matrix
#'
#' @return dgCMatrix.
#' @export


# Define a generic function for type conversion
setGeneric("convertToDgCMatrix", function(x) standardGeneric("convertToDgCMatrix"))

# Method for dgTMatrix
setMethod("convertToDgCMatrix", signature(x = "dgTMatrix"), function(x) {
  as(x, "dgCMatrix")
})

# Method for dense matrix
setMethod("convertToDgCMatrix", signature(x = "matrix"), function(x) {
  as(Matrix(x, sparse = TRUE), "dgCMatrix")
})