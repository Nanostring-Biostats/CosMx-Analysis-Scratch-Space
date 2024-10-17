
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
#' @export
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




#' for each cell, get the sum of x's values over its neighbors:
#' @param x A numeric vector
#' @param neighbors A (probably sparse) adjacency matrix
#' @importFrom Matrix rowSums
#' @export
neighbor_sum <- function(x, neighbors) {
  #Matrix::rowSums(t(t(neighbors != 0) * x))
  Matrix::rowSums(Matrix::t(Matrix::t(1*(neighbors != 0)) * x))
}


#' for each cell, get the mean of x's values over its neighbors:
#' @param x A numeric vector
#' @param neighbors A (probably sparse) adjacency matrix
#' @importFrom Matrix rowSums
#' @export
neighbor_mean <- function(x, neighbors) {
  neighbor_sum(x, neighbors) / pmax(Matrix::rowSums(neighbors != 0), 1)
}


#' for each cell, get the colSums of x over its neighbors:
#' @param x A matrix
#' @param neighbors A (probably sparse) adjacency matrix
#' @export
neighbor_colSums <- function(x, neighbors) {
  neighbors@x <- rep(1, length(neighbors@x))
  #neighbors <- Matrix::Diagonal(x=rep(1, nrow(neighbors)),names=rownames(neighbors)) %*% neighbors
  neighbors <- Matrix::Diagonal(x=rep(1, nrow(neighbors))) %*% neighbors
  neighbors@x[neighbors@x==0] <- 1
  out <- neighbors %*% x
  return(out)
}

#' for each cell, get the colMeans of x over its neighbors:
#' @param x A matrix
#' @param neighbors A (probably sparse) adjacency matrix
#' @export
neighbor_colMeans <- function(x, neighbors) {
  neighbors@x <- rep(1, length(neighbors@x))
  #neighbors <- Matrix::Diagonal(x=1/Matrix::rowSums(neighbors),names=rownames(neighbors)) %*% neighbors
  neighbors <- Matrix::Diagonal(x=1/Matrix::rowSums(neighbors)) %*% neighbors
  neighbors@x[neighbors@x==0] <- 1
  out <- neighbors %*% x
  return(out)
}

#' for each cell, tabulate the distinct values of x over its neighbors:
#' @param x A vector of categorical values
#' @param neighbors A (probably sparse) adjacency matrix
#' @export
neighbor_tabulate <- function(x, neighbors) {
  uniquevals <- unique(x)
  sapply(uniquevals, function(val){
    vec <- (x == val) * 1
    neighbor_sum(vec, neighbors)
  })
}

#' for each cell, tabulate the distinct values of x over its neighbors:
#' @param x A vector of categorical values
#' @param neighbors A (probably sparse) adjacency matrix
#' @export
neighbor_meantabulate <- function(x, neighbors) {
  uniquevals <- unique(x)
  sapply(uniquevals, function(val){
    vec <- (x == val) * 1
    neighbor_mean(vec, neighbors)
  })
}

#' Summarize a data frame's values over a neighborhood
#'
#' Numeric and character columns are respectively summed or tabulated across each cell's neighbors.
#' @param df A data frame, with rows aligned to the rows of neighbors. Can contain numeric and/or character columns. 
#' @param neighbors A neighbors adjacency matrix
#' @return A list with elements corresponding to data frame variables. Entries are either numeric vectors
#'  summarizing numeric variablies in conditionon or count matrices containing the tabulated values from character vectors.
#' @export
dataframe_neighborhood_summary <- function(df, neighbors) {
  out <- vector(mode = "list", length = ncol(df))
  for (i in 1:ncol(df)) {
    if (is.numeric(df[, i])) {
      out[[i]] <- neighbor_sum(df[, i], neighbors)
    }
    if (is.factor(df[, i])) {
      df[, i] <- as.character(df[, i])
    }
    if (is.character(df[, i])) {
      if (length(unique(df[, i])) > 200) {
        warning(paste0("Character variable ", colnames(df)[i], " has over 200 unique values. It is actually a numeric variable?"))
      }
      out[[i]] <- neighbor_tabulate(df[, i], neighbors)
    }
  }
  names(out) <- colnames(df)
  return(out)
}

#' Calculate neighborhood expression
#'
#' Calculates the expression profile of each cell's neighborhood
#' @param counts Single cell expression matrix
#' @param neighbors A neighbors adjacency matrix
#' @return A matrix in the same dimensions as \code{counts}, giving the expression profile of each cell's neighborhood.
#' @export
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



#' Subsample each row of a nearest neighbors matrix by 50%, 
#'  assuming each row has exactly the same number of non-zero entries
#' @param neighbors sparse matrix of neighbor relationships
#' @param p Subsampling rate, with 1 meaning all neighbors are kept. 
#' @return The same neighbors matrix, with p proportion of the entries in each row set to 0. 
#' @importFrom Matrix rowSums
#' @export
subsampleNeighborsByRow <- function(neighbors, p = 0.5) {
  # checks:
  if (length(unique(Matrix::rowSums(neighbors != 0))) > 1) {
    stop("This function requires that every row of neighbors has the same number of non-zero entries.")
  }
  
  # random vector of which neighbors to keep vs. reset to 0 (want to sample the same number from each row):
  K <- sum(neighbors[1, ] != 0)
  temp <- rep(seq_len(nrow(neighbors)), each = K) + rnorm(K * nrow(neighbors), mean = 0, sd = 1e-5)
  o <- order(temp) %% K + 1
  keep <- o <= round(K * p)
  
  # reset non-kept entries to 0
  neighbors <- Matrix::t(neighbors)
  neighbors@x <- neighbors@x * keep
  neighbors <- Matrix::t(neighbors)
  
  return(neighbors)
}