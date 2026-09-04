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
