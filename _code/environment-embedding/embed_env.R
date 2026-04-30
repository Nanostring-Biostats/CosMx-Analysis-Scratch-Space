## Concept:
# The goal here is to fit and project embeddings of single cells' "cellular neighborhoods", 
# i.e. the details of the tissue landscape surrounding them. 
# To achieve this, we begin with a single cell embedding matrix encoding single cells' expression states. 
# We also define a series of sparse matrices "neighbors5" "neighbors50" etc encoding each cell's K nearest neighbors in xy space.
# These neighbors define a cell's "neighborhood" at different distances / K values.
# We summarize a neighborhood using random projection sketching; that is, we apply a randomly defined linear transformation of M columns
#  to the matrix of single cell embeddings of neighborhood cells, then report the mean and mean squared values of the transformed columns. 
#  This defines an embedding over a cellular neighborhood.
# To fully describe a cell's spatial context, we use expanding neighborhoods: first its 5 nearest neighbors, then 50, then 200.
# These larger-radius neighborhoods are embedded by random projection sketching of the previous level's neighborhood embedding.
# We end with the cbind() of the k=5, k=50 and k=200 neighborhoods. 


#' Fit single-cell embeddings by calling scPearsonPCA
#' @param counts Counts matrix, cells x genes
#' @param tot Vector of cells' total counts
fitSingleCellEmbedding <- function(counts, tot) {
  
  genefreq <- scPearsonPCA::gene_frequency(Matrix::t(counts)) ## gene frequency (across all cells)
  
  pcaobj <- scPearsonPCA::sparse_quasipoisson_pca_seurat(
    x = Matrix::t(counts),
    totalcounts = tot,
    grate = genefreq,
    scale.max = 10, ## PC's reflect clipping pearson residuals > 10 SDs above the mean pearson residual
    do.scale = TRUE, ## PC's reflect as if pearson residuals for each gene were scaled to have standard deviation=1
    do.center = TRUE ## PC's reflect as if pearson residuals for each gene were centered to have mean=0
  )
  
  # return grate values and PC weights:
}


#' Apply single-cell embeddings
#' @param counts Counts matrix, cells x genes
#' @param tot Vector of cells' total counts
#' @param grate Vector of per-gene expression values
#' @param wts Weights from a scPearsonPCA fit
#' @return Embedding of single cell expression profiles
embedSingleCells <- function(counts, tot, grate, wts) {
  
}


#' Fit an embedding of cellular neighborhoods
#' Compute an embedding capturing information about cells' neighborhoods. 
#'  At escalating neighborhood sizes, used random projection sketching to summarize neighborhood content. 
#'  At default settings, summarizes single cell embeddings over k=5 neighborhoods.
#'  Then summarizes k=5 embeddings over k=50 neighborhoods, and summarizes k=50 embeddings at the k=200 level.
#'  Each cell's final embedding is the union of its k=5, k=50 and k=200 embeddings. 
#' @param scembed Matrix of single cell embeddings
#' @param xy 2-column matrix of cells' positions
#' @param ks Neighborhood sizes. 
#' @param sketchdims Number of random projection sketching variables to produce in each embedding
#' @param name description
fitNeighborhoodEmbedding <- function(scembed, xy, ks = c(5, 50, 200), sketchdims = c(20, 20, 20)) {
  
  
}


#' Apply a previously-defined neighborhood embedding
#' @param scembed Matrix of single cell embeddings
#' @param xy 2-column matrix of cells' positions
#' @param fit Neighborhood embedding params as defined by fitNeighborhoodEmbedding()
embedNeighborhoods <- function(scembed, xy, fit) {
  
}



#' Get random projection sketch of a matrix
#' @param mat A matrix
#' @param neighbors A sparse matrix encoding nearest spatial neighbors, corresponding to the rows of mat
#' @param nfeats How many features to produce
#' @return A matrix of random projection sketch values
randomProjectionSketch <- function(mat, nfeats) {
  
}