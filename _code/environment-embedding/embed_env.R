
#' Obtain single-cell embeddings by calling scPearsonPCA
embedSingleCells <- function() {
  
}


#' Embed cellular neighborhood information
#' Compute an embedding capturing information about cells' neighborhoods. 
#'  At escalating neighborhood sizes, used random projection sketching to summarize neighborhood content. 
#'  At default settings, summarizes single cell embeddings over k=5 neighborhoods.
#'  Then summarizes k=5 embeddings over k=50 neighborhoods, and summarizes k=50 embeddings at the k=200 level.
#'  Each cell's final embedding is the union of its k=5, k=50 and k=200 embeddings. 
#' @param scembed Matrix of single cell embeddings
#' @param xy 2-column matrix of cells' positions
#' @param ks Neighborhood sizes. 
#' @param sketchdims Number of random projection sketching variables to produce in each embedding
embedNeighborhoods <- function(scembed, xy, ks = c(5, 50, 200), sketchdims = c(20, 20, 20)) {
  
  
  
}




#' Get random projection sketch of a matrix
#' @param mat A matrix
#' @param nfeats How many features to produce
randomProjectionSketch <- function(mat, nfeats) {
  
}