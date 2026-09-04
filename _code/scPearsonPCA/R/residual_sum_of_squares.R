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


#' Compute each cell's sum of squares after projecting Pearson residuals onto the first PCs
#' 
#' @description
#' 
#' Given a PCA result from `sparse_quasipoisson_pca_seurat()`, computes for each cell the squared norm
#' of the orthogonal projection of its (centered/scaled/clipped) Pearson residual vector onto the first
#' `npcs` principal component directions. Equivalent to `rowSums(cell.embeddings[, 1:npcs]^2)`, but
#' computed directly WITHOUT ever forming the (dense) Pearson residuals matrix or the full cell embeddings.
#'
#' @param x A sparse genes x cells counts matrix (same genes used to fit `pcaobj`, or a superset)
#' @param pcaobj a PCA result list, as returned by `sparse_quasipoisson_pca_seurat()`
#' @param npcs the number of (leading) PCs to project onto. Defaults to all PCs available in `pcaobj`.
#' @param totalcounts  The total UMI counts across all gene targets in each cell.  This can be pre-computed by the user and passed in cases where the counts matrix `x` contains only a 
#' subset (of, say, highly variable) genes to be used in PCA. 
#'
#' @return A named numeric vector (one value per cell) of summed squared PC scores across the first `npcs` PCs.
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
  
  rss <- Matrix::colSums(as.matrix(scores)^2)
  names(rss) <- colnames(x)
  return(rss) 
}
