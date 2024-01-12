#' @title Update cell typing results with spatial context or other alternative data
#' 
#' @description
#' Takes cell typing results, then updates it based on alternative data types, 
#' e.g. spatial context, morphology, or protein expression. Existing cell typing results are 
#' put into Insitutype's likelihood framework, which then can use alternative data
#' as a prior to be updated by the expression data to get a new posterior probability 
#' of cell type.
#' Performs this operation by 
#' \enumerate{
#' \item deriving cell type profiles using InSituType:::Estep(), 
#' \item assigning cells to "cohorts" (clusters) derived from their alternative data
#' \item  Inputing the output of steps (1) and (2) into InSituType::insitutype() to 
#'  re-calculate cell type. 
#' }
#' Paths for using alternative data:
#' \enumerate{
#' \item Input your own \code{cohort} vector
#' \item Input a matrix of alternative data (\code{altdata}) to be automatically clustered into cohorts
#' \item Input \code{xy} positions (and possibly \code{tissue}). Then cells will be clustered 
#'  into cohorts based on the expression pattern of their 50 nearest neighoring cells.
#' }
#' @param celltype Vector of cell type assignments to be updated
#' @param counts Counts matrix (or dgCMatrix), cells * genes.
#' @param neg Vector of mean negprobe counts per cell
#' @param cohort Vector of cells' cohort memberships. Output of a spatial clustering algorithm makes for good cohorts. 
#' @param altdata Matrix of cells' alternative data values
#' @param xy 2-column matrix of cells' xy positions. 
#' @param tissue Vector giving cells' tissue IDs. Used to separate tissue with overlapping xy coordinates.
#' @param nb_size The size parameter to assume for the NB distribution.
#' @param assay.type A string specifying which assay values to use.
#' @importFrom InSituType Estep insitutype fastCohorting
#' @importFrom irlba irlba
spatialUpdate <- function(celltype, counts, neg, 
                          cohort = NULL, altdata = NULL, xy = NULL, tissue = NULL,
                          nb_size = 10, assay.type = "rna") {
  ## check alternative data args:
  if ((is.null(cohort) * is.null(altdata)) & is.null(xy)) {
    stop("Must supply cohort, altdata or xy")
  }

  ## process alternative data, obtaining cohort vector:
  if (is.null(cohort)) {
    if (!is.null(altdata)) {
      # make altdata from cells' neighborhoods:
      if (is.null(tissue)) {
        tissue = 1
      }
      neighbors <- nearestNeighborGraph(x = xy[, 1], y = xy[, 2], N = 50, subset = tissue) 
      neighborexpression <- get_neighborhood_expression(counts = counts, neighbors = neighbors) 
      # save top 20 PCs:
      temp <- irlba::irlba(neighborexpression, nv = 20)
      altdata <- temp$u %*% diag(temp$d)
    }
    # cluster altdata to get cohort:
    cohort <- InSituType::fastCohorting(mat = altdata, 
                                        n_cohorts = NULL, 
                                        gaussian_transform = TRUE) 
  }
  
  ## derive reference profiles from initial cell type vector:
  profiles <- InSituType:::Estep(counts = counts, 
                                 clust = celltype, 
                                 neg = neg)
 
  ## Run supervised cell typing with InSituType
  res <- InSituType::insitutype(x = counts, 
                                neg = neg, 
                                reference_profiles = profiles,
                                update_reference_profiles = FALSE)
  
  return(res)
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

#' Calculate neighborhood expression
#'
#' Calculates the expression profile of each cell's neighborhood
#' @param counts Single cell expression matrix
#' @param neighbors A neighbors adjacency matrix
#' @return A matrix in the same dimensions as \code{counts}, giving the expression profile of each cell's neighborhood.
get_neighborhood_expression <- function(counts, neighbors) {
  
  # check:
  if (nrow(counts) != ncol(neighbors)) {
    stop("misalignment between nrow(counts) and ncol(neighbors)")
  }
  # get clust-specific environment expression
  env <- neighbor_colMeans(counts, neighbors)
  rownames(env) <- rownames(neighbors)
  env <- as.matrix(env)
  return(env)
}

#' for each cell, get the colMeans of x over its neighbors:
#' @param x A matrix
#' @param neighbors A (probably sparse) adjacency matrix
neighbor_colMeans <- function(x, neighbors) {
  neighbors@x <- rep(1, length(neighbors@x))
  #neighbors <- Matrix::Diagonal(x=1/Matrix::rowSums(neighbors),names=rownames(neighbors)) %*% neighbors
  neighbors <- Matrix::Diagonal(x=1/Matrix::rowSums(neighbors)) %*% neighbors
  neighbors@x[neighbors@x==0] <- 1
  out <- neighbors %*% x
  return(out)
}

#' for each cell, get the colSums of x over its neighbors:
#' @param x A matrix
#' @param neighbors A (probably sparse) adjacency matrix
neighbor_colSums <- function(x, neighbors) {
  neighbors@x <- rep(1, length(neighbors@x))
  #neighbors <- Matrix::Diagonal(x=rep(1, nrow(neighbors)),names=rownames(neighbors)) %*% neighbors
  neighbors <- Matrix::Diagonal(x=rep(1, nrow(neighbors))) %*% neighbors
  neighbors@x[neighbors@x==0] <- 1
  out <- neighbors %*% x
  return(out)
}


