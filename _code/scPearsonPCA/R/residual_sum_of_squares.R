#' Compute each cell's sum of squared quasi-poisson Pearson residuals
#' 
#' @description
#' 
#' For each cell, computes `sum_g r_gc^2`, the sum over genes of squared Pearson residuals,
#' WITHOUT ever forming the (dense) Pearson residuals matrix. Useful as a per-cell diagnostic,
#' analogous to a goodness-of-fit statistic under the quasi-poisson null model.
#'
#' @param x A sparse genes x cells counts matrix 
#' @param totalcounts  The total UMI counts across all gene targets in each cell.  This can be pre-computed by the user and passed in cases where the counts matrix `x` contains only a 
#' subset (of, say, highly variable) genes to be used in PCA. 
#' @param grate The gene frequencies.  This can be pre-computed by the user and passed in cases where the counts matrix `x` contains only a 
#' subset (of, say, highly variable) genes to be used in PCA. 
#' @param quasi_poisson_variance_inflation factor of overdispersion for Variance relative to the mean.  i.e., V(x) = E(x) * phi, where phi is the variance inflation factor
#'
#' @return A named numeric vector (one value per cell) of summed squared Pearson residuals.
#'
#' @export  
residual_sum_of_squares <- function(x
                                    ,totalcounts = NULL
                                    ,grate = NULL
                                    ,quasi_poisson_variance_inflation = 1.01
){
  
  x <- check_x_is_dgcmatrix(x)
  phi <- check_phi(x, quasi_poisson_variance_inflation)
  totalcounts <- check_totalcounts(x, totalcounts)
  grate <- check_grate(x, grate, totalcounts)
  
  ### diagonal matrix (1/(phi*estimated gene frequency))
  root_grate_phi_diag <- Matrix::Diagonal(x = sqrt(1/(grate * phi)))
  dimnames(root_grate_phi_diag) <- list(names(grate), names(grate))
  
  ### diagonal matrix (1/totalcounts)
  root_tc_diag <- Matrix::Diagonal(x = sqrt(1/totalcounts))
  dimnames(root_tc_diag) <- list(colnames(x), colnames(x))
  
  ### y / sqrt(v(y)) 
  ytilde <- root_grate_phi_diag %*% x %*% root_tc_diag ## sparse genes x cells, same nnz as x
  
  root_grate_phi_row <- Matrix::Matrix(data = sqrt(grate / phi)
                                       ,nrow = 1
                                       ,dimnames = list(c(), names(grate)))
  
  ### r_gc = ytilde_gc - a_g*b_c, so sum_g(r_gc^2) = sum_g(ytilde_gc^2) - 2*b_c*sum_g(a_g*ytilde_gc) + b_c^2*sum_g(a_g^2)
  ytilde_sq <- ytilde
  ytilde_sq@x <- ytilde_sq@x^2
  
  term1 <- Matrix::colSums(ytilde_sq)
  term2 <- as.vector(root_grate_phi_row %*% ytilde)
  term3 <- sum(grate / phi)
  
  rss <- term1 - 2 * sqrt(totalcounts) * term2 + term3 * totalcounts
  names(rss) <- colnames(x)
  return(rss) 
}


#' Compute each cell's sum of squares in, and orthogonal to, the space of the first PCs
#' 
#' @description
#' 
#' Given a PCA result from `sparse_quasipoisson_pca_seurat()`, computes for each cell the squared norm
#' of the orthogonal projection of its (centered/scaled/clipped) Pearson residual vector onto the first
#' `npcs` principal component directions ("in-space"), as well as the squared norm of what's left over
#' after removing that projection ("orthogonal"). The two add up to the total sum of squares of the
#' (centered/scaled/clipped) Pearson residual vector. The in-space value is equivalent to
#' `rowSums(cell.embeddings[, 1:npcs]^2)`, but both values are computed directly WITHOUT ever forming
#' the (dense) Pearson residuals matrix or the full cell embeddings.
#'
#' @param x A sparse genes x cells counts matrix (same genes used to fit `pcaobj`, or a superset)
#' @param pcaobj a PCA result list, as returned by `sparse_quasipoisson_pca_seurat()`
#' @param npcs the number of (leading) PCs to project onto. Defaults to all PCs available in `pcaobj`.
#' @param totalcounts  The total UMI counts across all gene targets in each cell.  This can be pre-computed by the user and passed in cases where the counts matrix `x` contains only a 
#' subset (of, say, highly variable) genes to be used in PCA. 
#'
#' @return A list of 2 named numeric vectors (one value per cell each): `projection_sum_of_squares`
#' (summed squared PC scores across the first `npcs` PCs) and `orthogonal_sum_of_squares` (the
#' remaining sum of squares outside the span of those PCs).
#'
#' @export  
residual_projection_sum_of_squares <- function(x
                                               ,pcaobj
                                               ,npcs = NULL
                                               ,totalcounts = NULL
){
  
  feature.loadings <- 
    if(!is.null(pcaobj$feature.loadings)) pcaobj$feature.loadings else pcaobj$reduction.data@feature.loadings
  if(is.null(npcs)) npcs <- ncol(feature.loadings)
  stopifnot("'npcs' cannot exceed the number of PCs available in 'pcaobj'" = npcs <= ncol(feature.loadings))
  feature.loadings <- feature.loadings[, seq_len(npcs), drop = FALSE]
  
  x <- check_x_is_dgcmatrix(x)
  genes <- names(pcaobj$grate)
  stopifnot("'pcaobj$grate' must be a named vector" = !is.null(genes))
  stopifnot("genes used to fit 'pcaobj' are missing from 'x'" = all(genes %in% rownames(x)))
  x <- x[genes, , drop = FALSE]
  feature.loadings <- feature.loadings[genes, , drop = FALSE]
  
  phi <- check_phi(x, pcaobj$phi)
  totalcounts <- check_totalcounts(x, totalcounts)
  grate <- check_grate(x, pcaobj$grate, totalcounts)
  mean.pearson_residual <- pcaobj$mean.pearson_residual[genes]
  sd.pearson_residual <- pcaobj$sd.pearson_residual[genes]
  scale.max <- pcaobj$scale.max
  do.center <- pcaobj$do.center
  do.scale <- pcaobj$do.scale
  
  ### diagonal matrix (1/(phi*estimated gene frequency))
  root_grate_phi_diag <- Matrix::Diagonal(x = sqrt(1/(grate * phi)))
  dimnames(root_grate_phi_diag) <- list(genes, genes)
  
  ### diagonal matrix (1/totalcounts)
  root_tc_diag <- Matrix::Diagonal(x = sqrt(1/totalcounts))
  dimnames(root_tc_diag) <- list(colnames(x), colnames(x))
  
  ### y / sqrt(v(y)) 
  ytilde <- root_grate_phi_diag %*% x %*% root_tc_diag ## sparse genes x cells
  
  #### re-apply the same scale.max clipping used when 'pcaobj' was fit, so scores match the fitted PCs
  if(scale.max < Inf){
    clip_vals <- scale.max * sd.pearson_residual + mean.pearson_residual
    xvec <- ytilde@x   
    names(xvec) <- rownames(ytilde)[(ytilde@i + 1)] 
    cellnames_max_val_vec <- colnames(ytilde)[rep.int(seq_len(ncol(ytilde)), diff(ytilde@p))]
    max_val_vec <- clip_vals[names(xvec)] + 
      sqrt(grate[names(xvec)]/phi) * sqrt(totalcounts[cellnames_max_val_vec])
    ytilde@x <- pmin(xvec, max_val_vec) 
    rm(list = c("xvec", "max_val_vec", "cellnames_max_val_vec"))
  }
  
  ### fold the per-gene 1/sd scaling into the loadings, so only a npcs x genes by genes x cells product is ever formed
  root_loadings <- if(do.scale) feature.loadings / sd.pearson_residual else feature.loadings
  root_grate_phi_row <- Matrix::Matrix(data = sqrt(grate / phi), nrow = 1, dimnames = list(c(), genes))
  
  ### npcs x cells: PC scores computed directly, without ever forming the genes x cells residual matrix
  scores <- Matrix::crossprod(root_loadings, ytilde) - 
    as.vector(root_grate_phi_row %*% root_loadings) %o% sqrt(totalcounts)
  if(do.center){
    scores <- scores - as.vector(Matrix::crossprod(root_loadings, mean.pearson_residual))
  }
  
  projection_rss <- Matrix::colSums(as.matrix(scores)^2)
  names(projection_rss) <- colnames(x)
  
  ### total sum of squares of the SAME (centered/scaled/clipped) residual r'_gc = U_gc - A_g*b_c - M_g,
  ### using the identical building blocks as 'scores' above but without projecting onto 'root_loadings' first
  if(do.scale){
    inv_sd_diag <- Matrix::Diagonal(x = 1 / sd.pearson_residual)
    dimnames(inv_sd_diag) <- list(genes, genes)
    U <- inv_sd_diag %*% ytilde ## sparse, same nnz as ytilde
    A <- sqrt(grate / phi) / sd.pearson_residual
    M <- if(do.center) mean.pearson_residual / sd.pearson_residual else rep(0, length(genes))
  } else {
    U <- ytilde
    A <- sqrt(grate / phi)
    M <- if(do.center) mean.pearson_residual else rep(0, length(genes))
  }
  names(A) <- names(M) <- genes
  A_row <- Matrix::Matrix(data = A, nrow = 1, dimnames = list(c(), genes))
  M_row <- Matrix::Matrix(data = M, nrow = 1, dimnames = list(c(), genes))
  
  U_sq <- U
  U_sq@x <- U_sq@x^2
  
  b <- sqrt(totalcounts)
  term1 <- Matrix::colSums(U_sq)
  term2 <- as.vector(A_row %*% U)
  term3 <- as.vector(M_row %*% U)
  term4 <- sum(A^2)
  term5 <- sum(A * M)
  term6 <- sum(M^2)
  
  total_rss <- term1 - 2 * b * term2 - 2 * term3 + b^2 * term4 + 2 * b * term5 + term6
  names(total_rss) <- colnames(x)
  
  ### valid by Pythagoras since 'feature.loadings' columns are orthonormal (they come from RSpectra::svds)
  orthogonal_rss <- total_rss - projection_rss
  
  return(list(projection_sum_of_squares = projection_rss
             ,orthogonal_sum_of_squares = orthogonal_rss)) 
}
