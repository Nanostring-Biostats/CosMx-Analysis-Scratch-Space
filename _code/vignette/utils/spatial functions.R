# Copyright ©2024 Bruker Spatial Biology, Inc. All rights reserved. Subject to additional license terms and conditions provided separately by Bruker Spatial Biology, Inc.

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


#' for each cell, get the sum of x's values over its neighbors:
#' @param x A numeric vector
#' @param neighbors A (probably sparse) adjacency matrix
#' @importFrom Matrix rowSums
neighbor_sum <- function(x, neighbors) {
  #Matrix::rowSums(t(t(neighbors != 0) * x))
  Matrix::rowSums(Matrix::t(Matrix::t(1*(neighbors != 0)) * x))
}


#' for each cell, get the mean of x's values over its neighbors:
#' @param x A numeric vector
#' @param neighbors A (probably sparse) adjacency matrix
#' @importFrom Matrix rowSums
neighbor_mean <- function(x, neighbors) {
  neighbor_sum(x, neighbors) / pmax(Matrix::rowSums(neighbors != 0), 1)
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

