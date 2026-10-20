#' Compute gene signature scores as linear combinations of quasi-poisson Pearson residuals
#' 
#' @description
#' 
#' For a set of gene signatures (weighted linear combinations of genes, given by `wts`), computes the
#' cell-level score `sum_g w_g * r_gc`, where `r_gc` is the quasi-poisson Pearson residual of gene g in
#' cell c, WITHOUT ever forming the (dense) Pearson residuals matrix.
#'
#' @param x A sparse genes x cells counts matrix 
#' @param wts A genes x signatures matrix of weights, with rownames giving the gene names used in each signature 
#' @param totalcounts  The total UMI counts across all gene targets in each cell.  This can be pre-computed by the user and passed in cases where the counts matrix `x` contains only a 
#' subset (of, say, highly variable) genes to be used in PCA. 
#' @param grate The gene frequencies.  This can be pre-computed by the user and passed in cases where the counts matrix `x` contains only a 
#' subset (of, say, highly variable) genes to be used in PCA. 
#' @param quasi_poisson_variance_inflation factor of overdispersion for Variance relative to the mean.  i.e., V(x) = E(x) * phi, where phi is the variance inflation factor
#'
#' @return A cells x signatures matrix of signature scores.
#'
#' @export  
pearson_signatures <- function(x
                               ,wts
                               ,totalcounts = NULL
                               ,grate = NULL
                               ,quasi_poisson_variance_inflation = 1.01
){
  
  stopifnot("'wts' must have rownames giving gene names" = !is.null(rownames(wts)))
  genes <- rownames(wts)
  missing_genes <- setdiff(genes, rownames(x))
  if(length(missing_genes) > 0){
    stop(paste0("The following genes in 'wts' are missing from 'x': "
               ,paste0(missing_genes, collapse = ", ")))
  }
  if(is.null(colnames(wts))) colnames(wts) <- paste0("signature", seq_len(ncol(wts)))
  
  x <- check_x_is_dgcmatrix(x)
  x <- x[genes, , drop = FALSE]
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
  ytilde <- root_grate_phi_diag %*% x %*% root_tc_diag ## sparse genes x cells
  
  ### signatures x cells: weighted sum of y/sqrt(v(y)) across genes, for every signature at once
  part_a <- Matrix::crossprod(wts, ytilde) 
  
  ### signatures x cells: weighted sum of muhat/sqrt(v(y)) across genes, for every signature at once
  part_b <- as.vector(Matrix::crossprod(wts, sqrt(grate/phi))) %o% sqrt(totalcounts)
  
  embed <- as.matrix(Matrix::t(part_a - part_b))
  dimnames(embed) <- list(colnames(x), colnames(wts))
  return(embed) 
}
