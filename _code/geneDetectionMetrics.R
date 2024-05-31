

#' Wrapper function to call genes on all metrics:
#' @param counts Matrix of single cell gene expression profiles,  cells x genes
#' @param negcounts Matrix of single cell negative control expression profiles,  cells x negative control probes
#' @param xy 2-column matrix of cell xy positions
#' @return A vector of global p-values representing the ultimate word on a gene's utility, and 
#'         a matrix of detection p-values for genes x 4 metrics based on:
#' \itemize{
#'  \item Mean expression
#'  \item Extreme expression values
#'  \item Evidence of correlation with neighbors in physical space
#'  \item Evidence of correlation with neighbors in expression space 
#' }
#' @export
calcDetectionMetrics <- function(counts, negcounts, xy, subset_size = 100000) {
  
  # choose a subset:
  sub <- sample(seq_len(nrow(counts)), min(subset_size, nrow(counts)))
  
  # get normalized data:
  totcounts <- pmax(Matrix::rowSums(counts), 10) / 1000
  norm <- sweep(counts, 1, totcounts, "/")
  normneg <- sweep(negcounts, 1, totcounts, "/")
  
  # calc the metrics:
  ps <- matrix(NA, ncol(counts), 4,
               dimnames = c("mean", "extreme", "spatialcor", "expressioncor"))
  ps[, "mean"] <- calcMeanMetric()
  ps[, "extreme"] <- calcMeanMetric()
  ps[, "spatialcor"] <- calcSpatialCorMetric(norm = norm, normneg = normneg, xy = xy, sub = sub) 
  ps[, "expressioncor"] <- calcExpressionCorMetric(norm = norm, normneg = normneg, xy = xy, sub = sub)
  
  # get a global p-value using Fisher's method:
  teststats <- -2 * rowSums(log(ps))
  global_ps <- pchisq(teststats, df = 2 * ncol(ps), lower.tail = FALSE)
  
  out = list(global_ps, ps)
  return(out)
}



#' Score mean expression vs. what poisson would imply:
calcMeanMetric <- function(counts, negcounts) {
  
  meancounts <- Matrix::colMeans(counts)
  meannegcounts <- Matrix::colMeans(negcounts)
  
  ps <- 1 - pnorm(meancounts, mean = mean(meannegcounts), sd = sd(meannegcounts))
  return(ps)
}


#' Score extreme expression vs. what poisson would imply:
calcMaxMetric <- function(counts, negcounts) {
  
  # get expected negprobe rate:
  totcounts <- Matrix::rowSums(counts)
  meanegs <- Matrix(rowMeans(negcounts))
  negspertotcount <- sum(meanegs) / sum(totcounts)
  
  # get poission parameters for background in each cell:
  lambdas <- totcounts * negspertotcount
  
  # for each gene, compare observed counts vs. lambdas. Compute only for high expression values, and bonferroni correct:
  ps <- apply(counts, 2, function(x) {
    highcells <- (x > quantile(x, 0.999))
    tempps <- 1-ppois(x[highcells], lambda = lambdas[highcells]) 
    return(pmin(length(x) * tempps[order(tempps)[3]], 1)) # return the 3rd best p-value, after bonferroni correction
  })
}


#' Score correlation with spatial neighbors:
#' @param norm Normalized expression, cells * genes
#' @param normneg Negative control data normalized in the same manner as the genes in \code{norm}. 
#' @param xy 2-column matrix of cell xy positions
#' @param sub Vector giving the subset of cells to be considered
calcSpatialCorMetric <- function(norm, normneg, xy, sub = TRUE) {
  # get the each cell's 50-cell neighborhood:
  neighbors <- nearestNeighborGraph(x = xy[, 1], y = xy[, 2], N = k, sub = tissue)
  
  # get the mean profile over each neighborhood:
  neighbormean <- get_neighborhood_expression(counts = norm,
                                              neighbors = neighbors[sub, ])
  negneighbormean <- get_neighborhood_expression(counts = normneg,
                                              neighbors = neighbors[sub, ])
  
  # get the correlation between single cell and neighborhood expression for each gene:
  cors <- fastCor(norm, neighbormean)
  negcors <- fastCor(negnorm, negneighbormean)
  
  # compare the gene correlations to the negative control correlations to get a p-value:
  # (open question: use quantiles, or assume the neg values are normally distributed? The latter would be ideal if it was true. Test it.)
  
  return(ps)
}



#' Score correlation with expression space neighbors:
#' @param norm Normalized expression, cells * genes
#' @param normneg Negative control data normalized in the same manner as the genes in \code{norm}. 
#' @param xy 2-column matrix of cell xy positions
#' @param sub Vector giving the subset of cells to be considered
calcExpressionCorMetric <- function(norm, normneg, xy, sub = TRUE) {
  
  # project to 20 PCs:
  temp <- irlba::irlba(norm[sub, ], nv = 20)
  pcs <- temp$u %*% diag(temp$d)
  
  # get nearest neighbors:
  expressionneighbors <- FNN::get.knnx(data = pcs, query = pcs, k = 26)$nn.index
  
  # get mean expression in nearest PC neighbors:
  neighbormean <- get_neighborhood_expression(counts = norm,
                                              neighbors = expressionneighbors[, -1])
  negneighbormean <- get_neighborhood_expression(counts = normneg,
                                              neighbors = expressionneighbors[, -1])

  # get the correlation between single cell and neighborhood expression for each gene:
  cors <- fastCor(norm, neighbormean)
  negcors <- fastCor(negnorm, negneighbormean)
  
  # compare the gene correlations to the negative control correlations to get a p-value:
  # (open question: use quantiles, or assume the neg values are normally distributed? The latter would be ideal if it was true. Test it.)
  
  return(ps)
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


#' Efficient use of matrix algebra to get correlation between columns of 2 matrices:
fastCor <- function(A, B) {
  
  # center and scale:
  A_centered <- sweep(A, 2, Matrix::colMeans(A))
  B_centered <- sweep(B, 2, Matrix::colMeans(B))
  A_norm <- sweep(A_centered, 2, sqrt(colSums(A_centered^2)), "/")
  B_norm <- sweep(B_centered, 2, sqrt(colSums(B_centered^2)), "/")
  
  # Compute the dot product of normalized vectors
  correlations <- colSums(A_norm * B_norm)
  return(correlations)
}
# test:
if (FALSE) {
  mat1 = matrix(rnorm(1000),100)
  mat2 = matrix(rnorm(1000),100)
  expect_true(max(abs(fastCor(mat1, mat2) - diag(cor(mat1, mat2)))) < 1e-6)
  expect_true(cor(fastCor(mat1, mat2), diag(cor(mat1, mat2))) > 0.99)
}
