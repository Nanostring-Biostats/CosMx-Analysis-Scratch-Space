
#' calculate neighborhood expression
#' @param x matrix to get neighbor colMeans over
#' @param neighbors Sparse matrix giving neighbor relationships
#' @param makedense Logical, for whether to convert to dense matrix
#' @return A matrix in the dimension of x holding colMeans of x over neighbors
#' @export
getNeighborhoodExpression <- function(x, neighbors, makedense = TRUE) {
  # calculate colMeans of expression in the neighborhood:
  neighbors@x <- rep(1, length(neighbors@x))
  neighbors <- Matrix::Diagonal(x=1/Matrix::rowSums(neighbors)) %*% neighbors
  neighbors@x[neighbors@x==0] <- 1
  out <- neighbors %*% x
  if (makedense) {  #<------------- converting to dense matrix, since expectation is dense data
    out <- as.matrix(out)
  }
  return(out) 
}

#' function to get neighbors via all the various paths:
#' @param xy Cells' xy coordinates - a 2-column matrix aligned to the rows of mat
#' @param neighbors Sparse matrix of neighbor relationships. If not provided, xy will be used to calculate thi]s.
#' @param tissue Vector of cells' tissue IDs
#' @param controlnames Tissue IDs of control tissue(s)
#' @param k Number of neighbors to use to define a neighborhood
#' @param radius Radius used to define neighborhoods (choose this or K)
#' @param verbose TRUE to get more messages
#' @return sparse matrix holding neighbor assignments
#' @export
getNeighbors <- function(xy, neighbors, tissue, controlnames, k, radius, verbose) {
  # need xy or neighbors:
  if (is.null(xy) && is.null(neighbors)) {
    stop("need to provide either xy coords or a neighbors network")
  }
  # if xy AND neighbors provided, use neighbors and not xy:
  if (!is.null(xy) && !is.null(neighbors)) {
    message("both xy and neighbors were provided. Using the neighbors provided rather than recalculating.")
    xy <- NULL
  }
  # check the integrity of the neighbors object:
  if (!is.null(neighbors)) {
    if (!identical(dim(neighbors), rep(nrow(counts), 2))) {
      stop("neighbors object must be a matrix or sparse matrix with nrow and ncol both equal nrow(counts)")
    }
    propnonzero <- (sqrt(Matrix::nnzero(neighbors)) / (dim(neighbors)[1]))^2
    if (propnonzero > 0.05) {
      warning("The neighbors matrix provided isn't very sparse - this could indicate an error.")
    }
    
  }
  # if calculating neighbors, then need k or radius, and not both
  if (!is.null(xy)) {
    if (is.null(k) && is.null(radius)) {
      stop("must provide either k or radius for neighbors to be calculated")
    }
    if (!is.null(k) && !is.null(radius)) {
      message("both k and radius were provided. Proceeding with k; radius will have no impact.")
      radius <- NULL
    }
  }
  
  ## build neighbors graph:
  if (is.null(neighbors)) {
    if (verbose) {
      print("building nearest neighbors network")
    }
    if (is.null(tissue)) {
      tissue = 1
    }
    if (!is.null(k)) {
      neighbors <- nearestNeighborGraph(x = xy[, 1], y = xy[, 2], N = k, subset = tissue)
    }
    if (!is.null(radius)) {
      neighbors <- radiusBasedGraph(x = xy[, 1], y = xy[, 2], R = radius, subset = tissue)
    }
    rownames(neighbors) <- rownames(xy)
    colnames(neighbors) <- rownames(xy)
  }
  return(neighbors)
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

#' Create spatial network from neighbors within radius R
#'
#' For each cell identify neighbors within distance \code{R} in Euclidean space
#' and create an edge between them in graph structure, optionally subset cells
#' (see Details).
#'
#' Edges will only be created for cells that have the same \code{subset} value,
#' usually the slide column id but could also be a slide plus FOV id to only
#' create edges within an FOV.
#'
#' @param x spatial coordinate
#' @param y spatial coordinate
#' @param R radius
#' @param subset same length as x,y (see Details)
#'
#' @return sparse adjacency matrix with distances
#' @importFrom data.table data.table
#' @importFrom data.table rbindlist
#' @importFrom Matrix sparseMatrix
#' @importFrom spatstat.geom ppp
#' @importFrom spatstat.geom closepairs
radiusBasedGraph <- function(x, y, R, subset=1) {
  DT <- data.table::data.table(x = x, y = y, subset = subset)
  radiusNeighbor <- function(i) {
    subset_dt <- DT[subset == i]
    idx <- which(DT[["subset"]] == i)
    pp <- spatstat.geom::ppp(subset_dt$x, subset_dt$y,
                             range(subset_dt$x), range(subset_dt$y))
    cp <- spatstat.geom::closepairs(pp, R)
    ij <- data.table::data.table(i = idx[cp$i],
                                 j = idx[cp$j],
                                 x = cp$d)
    return(ij)
  }
  ij <- data.table::rbindlist(lapply(unique(subset), radiusNeighbor))
  adj.m <- Matrix::sparseMatrix(i = ij$i, j = ij$j, x = ij$x, dims = c(nrow(DT), nrow(DT)))
  return(adj.m)
}