message("Key functions:
      - findSafeGenes: identify genes safe from excessive bias from segmentation errors.
      - getSubclusteringGenes: identify highly variable genes in a counts matrix, for use in clustering.")


#' Identify genes useful for unsupervised clustering of closely-related cell types 
#' Intended for subclustering of a broad cell type, e.g. of immune cells into T-cells, B-cells, etc... 
#' @param mat Expression matrix from just the cell type in question, cells in rows and genes in columns.
#' @param xy Spatial coordinates of cells 
#' @param varratiothresh Genes with observed variance / predicted variance above this threshold will be called HVGs. Higher = a more stringent filter.
#' @param expressionthresh Only keep genes with average raw counts above this level in your cell type of interest. 
#' @param loess.span Parameter guiding loess fit used in hvg search.
#' @return A vector of gene named for use in unsupervised subclustering
getSubclusteringGenes <- function(mat, varratiothresh = 1, expressionthresh = 0.2, loess.span = 0.3) {
  
  ## identify hvgs in the counts matrix:
  RawTargetVar <- colvars(mat)
  RawTargetMean <- Matrix::colMeans(mat)
  use <- RawTargetVar > 0
  loessfit <- loess(log10(RawTargetVar[use]) ~ log10(RawTargetMean[use]),
                    span = loess.span)
  varratios <- RawTargetVar[use] / 10^loessfit[["fitted"]]
  hvgs <- names(which(varratios > varratiothresh))
  
  ## identify genes with decent expression in your data:
  means <- Matrix::colMeans(mat)
  decentexpressers <- names(which(means > expressionthresh))
  
  keepgenes <- intersect(hvgs, decentexpressers)
  return(keepgenes)
} 


#' Find genes safe from major contamination from cell segmentation errors
#' @param mat Counts matrix for the whole dataset, cells in rows and genes in columns.
#' @param xy Spatial coordinates of cells 
#' @param ismycelltype Logical vector, for whether cells belong to your cell type of interest
#' @param tissue Optional vector of cells' tissue IDs, used in case tissues overlap in xy space
#' @return A list with two elements:
#' \itemize{
#'  \item self2neighborratio: A vector of genes' ratios between the cell type in 
#'  question and neighborhoods of the cell type in question. High values indicate
#'   lower risk of bias from segmentation errors
#'  \item safegenes: A vector of gene names passing the filter
#' }
findSafeGenes <- function(counts, xy, ismycelltype, tissue = NULL, Nneighbors = 50, self_vs_neighbor_threshold = 1.75) {
  
  # get spatial neighbors:
  if (is.null(tissue)) {
    subset = 1
  } else {
    subset = tissue
  }
  neighbors <- nearestNeighborGraph(x = xy[, 1], y = xy[, 2], N = Nneighbors, subset=subset)

  # disqualify neighbors from your cell type:
  neighbors[, ismycelltype] <- 0
  
  # for just your cell type of interest, get total neighbor expression from other cell types
  totalenvcounts <- neighbor_colSums(x = counts, neighbors = neighbors[ismycelltype, ])
  meanenv <- Matrix::colMeans(totalenvcounts) 
  meanenv <- meanenv / Nneighbors
  # get mean expression 
  meanself <- Matrix::colMeans(counts[ismycelltype, ])
  
  out <- list(self2neighborratio = meanself / meanenv,
              safegenes = names(which((meanself / meanenv) > self_vs_neighbor_threshold)))
  return(out)
  
}



#' Create spatial network from N nearest neighbors
#'
#' For each cell identify \code{N} nearest neighbors in Euclidean space and
#' create an edge between them in graph structure, optionally subset cells (see
#' Details).
#'
#' Edges will only be created for cells that have the same \code{subset} value,
#' usually the slide column id but could also be a slide plus FOV id to only
#' create edges within an FOV.
#'
#' @param x spatial coordinate
#' @param y spatial coordinate
#' @param N number of nearest neighbors
#' @param subset same length as x,y (see Details)
#'
#' @return sparse adjacency matrix with distances
#' @importFrom data.table data.table
#' @importFrom data.table rbindlist
#' @importFrom spatstat.geom nnwhich
#' @importFrom spatstat.geom nndist
#' @importFrom Matrix sparseMatrix
#' @export
nearestNeighborGraph <- function(x, y, N, subset=1) {
  DT <- data.table::data.table(x = x, y = y, subset = subset)
  nearestNeighbor <- function(i) {
    subset_dt <- DT[subset == i]
    idx <- which(DT[["subset"]] == i)
    ndist <- spatstat.geom::nndist(subset_dt[, .(x, y)],
                                   k=1:N)
    nwhich <- spatstat.geom::nnwhich(subset_dt[, .(x, y)],
                                     k=1:N)
    ij <- data.table::data.table(i = idx[1:nrow(subset_dt)],
                                 j = idx[as.vector(nwhich)],
                                 x = as.vector(ndist))
    return(ij)
  }
  ij <- data.table::rbindlist(lapply(unique(subset), nearestNeighbor))
  adj.m <- Matrix::sparseMatrix(i = ij$i, j = ij$j, x = ij$x, dims = c(nrow(DT), nrow(DT)))
  return(adj.m)
}

#' for each cell, get the sum of x's values over its neighbors:
#' @param x A numeric vector
#' @param neighbors A (probably sparse) adjacency matrix
#' @importFrom Matrix rowSums
#' @export
neighbor_sum <- function(x, neighbors) {
  Matrix::rowSums(Matrix::t(Matrix::t(1*(neighbors != 0)) * x))
}

#' for each cell, get the colSums of x over its neighbors:
#' @param x A matrix
#' @param neighbors A (probably sparse) adjacency matrix
#' @export
neighbor_colSums <- function(x, neighbors) {
  neighbors@x <- rep(1, length(neighbors@x))
  neighbors <- Matrix::Diagonal(x=rep(1, nrow(neighbors))) %*% neighbors
  neighbors@x[neighbors@x==0] <- 1
  out <- neighbors %*% x
  return(out)
}


#' Calculate column standard deviations without converting to dense matrix.
#' @param x A sparse dgCMatrix
#' @examples 
#'  
#' set.seed(30)
#' sm <- Matrix::sparseMatrix(i=sample(1:30,450,replace=TRUE),
#'                            j=sample(1:20,450,replace=TRUE),x=runif(450))
#' all.equal(colSDs(sm), apply(sm, 2, sd))
#' 
#' @export
colvars <- function(x){
  
  ## std = \sqrt((\sum x_i^2 - n \bar{x}^2)/(n-1))
  if(inherits(x, "sparseMatrix")){
    x2 <- x
    x2@x <- x2@x^2
    x2.sums <- Matrix::colSums(x2)
    x.means <- Matrix::colMeans(x)
    num <- pmax((x2.sums - nrow(x)*x.means^2 ),0)  ## floating point error can cause this to be negative i.e, -4e-11, when sd is 0.
    colvars <- (num / (nrow(x) - 1))
  } else {
    colvars <- apply(x, 2, var)
  }
  
  return(colvars) 
}
  


