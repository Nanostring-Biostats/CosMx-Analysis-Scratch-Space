
#' @title moransMultiFast
#' @description Function for calculating Moran's I and Pvals for multiple genes using matrices
#'
#' Given a spatial weights matrix and gene expression matrix, returns Moran's Is and estimated Pvals
#' @param X gene expression matrix corresponding to selected cells and genes
#' @param W spatial weights matrix defining neighbors for selected cells
#' @param perm number of permutations used in MC estimation
#' 
#' @return morans_res with I, pval, and padj for each gene
#' @export
moransMultiFast <- function(X, W, perm, compute_sparse = TRUE, compute_perm = FALSE) {
  
  mean_of_cols <- sparseMatrixStats::colMeans2(X, useNames = T)
  sd_of_cols <- sparseMatrixStats::colSds(X, useNames = T)
  ncells <- nrow(X) 
  
  moransMatFast <- function(Z_index=NULL, Z=Z, m1C1=miC1, W=W) {
    # colSums(Z * (W %*% Z)) is equivalent to diag(t(Z) %*% W %*% Z) but less computationally intensive
    if(is.null(Z_index)){
      I <- Matrix::colSums(Z * (W %*% Z))/m1C1
    } else {
      Z_shuffle <- Z[Z_index,]
      I <- Matrix::colSums(Z_shuffle * (W %*% Z_shuffle))/m1C1 
      rm(Z_shuffle)
      gc()
    }
    return(I)
  }
  
  if(compute_sparse){
    M <- Matrix::sparseMatrix(i=1:nrow(W),j=1:nrow(W),x=1/Matrix::rowSums(W)
                              ,dimnames = dimnames(W))
    Wrs <- M %*% W
    WX <- Wrs %*% X 
    lag_mean_of_cols <- sparseMatrixStats::colMeans2(WX, useNames = T)
    diag_X_lagX <- Matrix::colSums(X * WX)
    Istat <- (diag_X_lagX - ncells * (lag_mean_of_cols * mean_of_cols)) / 
      # ((ncells-1)*sd_of_cols^2) ## ncells - 1 is technically correct, use ncells to match previous
      ((ncells)*sd_of_cols^2)
    
  } else {
    Z = Matrix::t( (Matrix::t(X) - mean_of_cols) / sd_of_cols )
    rownames(Z) <- rownames(X)
    colnames(Z) <- colnames(X)
    
    # Z_2 <- scale(as.matrix(X_dgC)) # z transform gene expression matrix on a per gene basis
    m1 <- rep(1,nrow(Z)) # Create vector of 1s with length equal to the number of cells
    m1C1 <- sum(t(m1) %*% W %*% m1) # calculate denominator for vectorized version of moran's I for multiple genes
    
    Istat <- moransMatFast(Z_index = NULL, Z, m1C1, W)   
  }
  
  if(compute_perm){
    # Calculate moran's I for all genes in expression matrix
    #spat_moran <- moransMatFast(Z_index = NULL, Z, m1C1, W) 
    spat_moran_sample_indices <- replicate(perm, sample(1:nrow(Z)))
    # Monte carlo simulation for n permutations of a randomized expression matrix
    num_threads <- getOption("MatrixExtra.nthreads")
    options("MatrixExtra.nthreads" = parallel::detectCores()-1)
    spat_moran_mc <- sapply(1:perm,function(i)moransMatFast(spat_moran_sample_indices[,i],Z,m1C1,W)) 
    options("MatrixExtra.nthreads" = num_threads)
    
    ## FOR LOOP SOLUTION FOR SAME PROCESS
    # spat_moran_mc_vec <- vector(length = perm*length(spat_moran))
    # for(i in 1:perm){
    #   temp_vec <- moransMatFast(spat_moran_sample_indices[,i],Z, m1C1, W)
    #   spat_moran_mc_vec[((i-1)*length(spat_moran)) + 1:length(spat_moran)] <- temp_vec
    #   rm(temp_vec)
    #   }
    # spat_moran_mc <- Matrix::Matrix(spat_moran_mc_vec, ncol=perm, byrow=T)
    
    EIsim <- rowSums(spat_moran_mc)/perm # Calculate average or estimate I values from MC simulation I values
    seIsim <- matrixStats::rowSds(spat_moran_mc) # Calculate std dev for simulated I values
    zsim <- (Istat - EIsim)/seIsim # Use estimated Is and simulated I std devs to z-transform actual I values
    pz <- 1 - pnorm(zsim) # calculate probability of simulated Is being greater than actual Is assuming simulated Is follow normal distributions
    
    spat_moran_res <- data.table::data.table(gene=names(Istat), I = Istat, pval = pz) 
    spat_moran_res[,padj:=p.adjust(pval, method="hochberg")]
    
  }
  if(!compute_perm){
    
    c1 <- Matrix::rowSums(Wrs)
    S0 <- sum(c1)
    S1 <- sum((Wrs * Wrs) + (Wrs * Matrix::t(Wrs)))
    S2 <- sum((Matrix::rowSums(Wrs) + Matrix::colSums(Wrs))^2)
    n <- nrow(X)
    consts <- 
      list(n = n, n1 = n-1, n2 = n-2, n3 = n-3, nn = n**2, S0 = S0, 
           S1 = S1, S2 = S2)
    
    K <- colKurtosis(X) 
    expected <- (-1)/consts$n1
    
    
    ### asymptotic version of variance
    S02 <- consts$S0 * consts$S0 
    
    VI <- consts$n * (consts$S1 * (consts$nn - 3 * consts$n + 3) - consts$n * 
                        consts$S2 + 3 * S02)
    tmp <- K * (consts$S1 * (consts$nn - consts$n) - 2 * consts$n * consts$S2 + 
                  6 * S02)
    if (any(tmp > VI)){
      msg <- paste0(names(tmp)[tmp > VI],collapse=",")
      msg <- paste0("Kurtosis overflow for genes: "
                    ,msg
                    ,"\ndistribution of variable does not meet test assumptions")
      warning(msg)
    }
    
    VI <- (VI - tmp)/(consts$n1 * consts$n2 * consts$n3 * S02)
    VI <- VI - expected^2
    ZI <- (Istat - expected) / sqrt(VI)
    
    spat_moran_res <- data.table::data.table(gene=names(ZI), I = Istat, pval = 1-pnorm(ZI), expected_I = expected, seI = sqrt(VI)) 
    spat_moran_res[,padj:=p.adjust(pval, method="hochberg")]
    
  }
  return(spat_moran_res) 
}
#' @title colKurtosis
#' @description kurtosis by column for sparse matrices
colKurtosis <- function(sm){
  n <- nrow(sm)
  sm4 <- sm3 <- sm
  sm4@x <- sm4@x ** 4
  sm3@x <- sm3@x ** 3
  ex <- Matrix::colMeans(sm)
  sigmas2 <- sparseMatrixStats::colSds(sm)^2 * (n-1) / n
  num <- Matrix::colMeans(sm4) - 4* ex * Matrix::colMeans(sm3) + 
    6*ex^2*sigmas2 + 3*ex**4
  denom <- sigmas2**2 
  return(num / denom)
}
