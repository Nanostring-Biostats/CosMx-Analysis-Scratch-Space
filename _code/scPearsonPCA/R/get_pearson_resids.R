#' Compute the full dense quasi-poisson Pearson residuals matrix
#' 
#' @description
#' 
#' Directly computes the (dense) genes x cells Pearson residuals matrix for counts matrix `x`,
#' without the genes x genes algebra trick used by `sparse_quasipoisson_pca_seurat()`.
#' Intended for validation on small matrices, since the resulting dense matrix can overwhelm
#' memory for large datasets.
#'
#' @param x A sparse genes x cells counts matrix, or a plain numeric vector of counts for a single gene (in which case `totalcounts` must be provided) 
#' @param totalcounts  The total UMI counts across all gene targets in each cell.  This can be pre-computed by the user and passed in cases where the counts matrix `x` contains only a 
#' subset (of, say, highly variable) genes to be used in PCA. 
#' @param grate The gene frequencies.  This can be pre-computed by the user and passed in cases where the counts matrix `x` contains only a 
#' subset (of, say, highly variable) genes to be used in PCA. 
#' @param quasi_poisson_variance_inflation factor of overdispersion for Variance relative to the mean.  i.e., V(x) = E(x) * phi, where phi is the variance inflation factor
#'
#' @return A dense genes x cells matrix of Pearson residuals.
#'
#' @export  
get_pearson_resids <- function(x
                               ,totalcounts = NULL
                               ,grate = NULL
                               ,quasi_poisson_variance_inflation = 1.01
){
  
  if(is.null(dim(x))){
    ### a single gene's own counts can't identify its gene frequency, so totalcounts is required here
    stopifnot("If 'x' is a plain vector (counts for a single gene), 'totalcounts' must be provided" = 
                !is.null(totalcounts))
    x <- Matrix::Matrix(matrix(x, nrow = 1, dimnames = list(NULL, names(x))), sparse = TRUE)
  }
  
  x <- check_x_is_dgcmatrix(x)
  if(ncol(x) <= 50 && is.null(totalcounts)){
    warning("Gene frequency is being computed from the observed counts in 'x' only. "
           ,"It will often make more sense to compute gene frequency from panel-wide total counts "
           ,"by providing a 'totalcounts' vector.")
  }
  phi <- check_phi(x, quasi_poisson_variance_inflation)
  totalcounts <- check_totalcounts(x, totalcounts)
  grate <- check_grate(x, grate, totalcounts)
  
  ### mu_gc = n_c * p_g, the quasi-poisson mean under the null model
  mu <- outer(grate, totalcounts)
  
  resids <- (as.matrix(x) - mu) / sqrt(phi * mu)
  dimnames(resids) <- dimnames(x)
  return(resids) 
}
