#' Compute the gene x gene correlation matrix of quasi-poisson Pearson residuals
#' 
#' @description
#' 
#' Computes the (dense) genes x genes Pearson correlation matrix of the quasi-poisson Pearson
#' residuals for counts matrix `x`, reusing the same genes x genes crossproduct algebra as
#' `sparse_quasipoisson_pca_seurat()`, WITHOUT ever forming the (dense) cells x genes residuals matrix.
#'
#' @param x A sparse genes x cells counts matrix 
#' @param totalcounts  The total UMI counts across all gene targets in each cell.  This can be pre-computed by the user and passed in cases where the counts matrix `x` contains only a 
#' subset (of, say, highly variable) genes to be used in PCA. 
#' @param grate The gene frequencies.  This can be pre-computed by the user and passed in cases where the counts matrix `x` contains only a 
#' subset (of, say, highly variable) genes to be used in PCA. 
#' @param quasi_poisson_variance_inflation factor of overdispersion for Variance relative to the mean.  i.e., V(x) = E(x) * phi, where phi is the variance inflation factor
#'
#' @return A dense genes x genes correlation matrix.
#'
#' @export  
pearson_cormat <- function(x
                          ,totalcounts = NULL
                          ,grate = NULL
                          ,quasi_poisson_variance_inflation = 1.01
){
  
  x <- check_x_is_dgcmatrix(x)
  nn <- ncol(x)
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
  
  ### gene-wise average of y/sqrt(v(y)) - muhat/sqrt(v(y))
  mean.pearson_residual <- Matrix::rowMeans(ytilde) - 
    mean(sqrt(totalcounts)) * sqrt(grate/phi)
  
  root_grate_phi_row <- Matrix::Matrix(data = sqrt(grate / phi)
                                       ,nrow = 1
                                       ,dimnames = list(c(), names(grate)))
  root_tc_col <- 
    Matrix::Matrix(sqrt(totalcounts), ncol = 1
                   ,dimnames = list(c(colnames(x)), c()))
  
  ytilde_muhat <- ytilde %*% root_tc_col ### g x 1
  ytilde_muhat <- ytilde_muhat %*% root_grate_phi_row ### g x g
  
  #### crossproduct of (uncentered) pearson residuals
  qp <- ytilde %*% Matrix::t(ytilde) - 
    ytilde_muhat - Matrix::t(ytilde_muhat) + 
    sum(root_tc_col^2) * Matrix::t(root_grate_phi_row) %*% root_grate_phi_row
  
  #### diag(qp) = sum_c(r_gc^2), giving the residual sd's the same way as sparse_quasipoisson_pca_seurat()
  sd.pearson_residual <- sqrt((Matrix::diag(qp) - nn * mean.pearson_residual^2) / (nn - 1))
  
  #### sum_c(r_gc*r_hc) - n*mean_g*mean_h = sum_c((r_gc-mean_g)(r_hc-mean_h)), i.e. (n-1)*Cov(g,h)
  qp_centered <- qp - nn * Matrix::tcrossprod(mean.pearson_residual)
  
  #### Cor(g,h) = Cov(g,h)/(sd_g*sd_h)
  inv_sd_diag <- Matrix::Diagonal(x = 1 / sd.pearson_residual)
  cormat <- as.matrix(inv_sd_diag %*% qp_centered %*% inv_sd_diag) / (nn - 1)
  dimnames(cormat) <- list(names(grate), names(grate))
  diag(cormat) <- 1 ## avoid floating-point drift off exactly 1
  return(cormat) 
}
