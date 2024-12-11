# dev notes:
# - contam filter needs ALL data; hvg filter needs just the selected cell type
# - implement both separately, w/o a wrapper. ppl can call both independently.
# - contam ratio function: 
# -- get neighbors
# -- redact neighbors of same cell type
# -- get env expression matrix
# -- get colMeans in self, in neighbors, take same ratio


#' Identify genes useful for supervised classification of closely-related cell types 
#' @param ref Reference matrix for the cell types in question
#' @param ratiothresh Only keep genes with at least this much of a ratio between the max and min cell types.
#' @param minquantilethresh Only keep genes above this quantile of the reference profile for at least one cell type. 
#' @return A vector of gene name for use in supervised cell subtyping
getSubtypingGenes <- function(ref, ratiothresh = 2, minquantilethresh = 0.5) {
  
  ## checks
  if (is.null(rownames(ref))) {
    stop("the reference profiles matrix (ref) needs row names")
  }
  
  ## identify hvgs in the ref profiles:
  mins <- apply(ref, 1, min)
  mins <- pmax(mins, min(mins[mins > 0], na.rm = TRUE))
  maxes <- apply(ref, 1, max)
  maxes <- pmax(maxes, min(mins))
  bigratios <- names(which((maxes / mins) > ratiothresh) )
  
  ## identify genes with decent expression in at least one of the refprofiles
  q <- sweep(ref, 2, apply(ref, 2, quantile, minquantilethresh), ">")
  decentexpressers <- rownames(q)[rowSums(q) > 0]
  
  keepgenes <- intersect(bigratios, decentexpressers)
  return(keepgenes)
} 




#' Identify genes useful for unsupervised clustering of closely-related cell types 
#' Intended for subclustering of a broad cell type, e.g. of immune cells into T-cells, B-cells, etc... 
#' @param mat Expression matrix from just the cell type in question, cells in rows and genes in columns.
#' @param xy Spatial coordinates of cells 
#' @param varratiothresh Genes with observed variance / predicted variance above this threshold will be called HVGs. Higher = a more stringent filter.
#' @param expressionthresh Only keep genes with average raw counts above this level in your cell type of interest. 
#' @return A vector of gene named for use in unsupervised subclustering
getSubclusteringGenes <- function(mat, loess.span = 0.3, varratiothresh = 1, expressionthresh = 0.2) {
  
  ## identify hvgs in the counts matrix:
  RawTargetVar <- apply(mat, 2, var)
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
findSafeGenes <- function(counts, xy, ismycelltype, tissue = NULL, Nneighbors = 50) {
  
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
              safegenes = names(which((meanself / meanenv) > 1)))
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
  #Matrix::rowSums(t(t(neighbors != 0) * x))
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


