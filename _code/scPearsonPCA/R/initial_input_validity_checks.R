

check_x_is_dgcmatrix <- function(x){
  if(!inherits(x, "sparseMatrix")){
    stop("'x' should be of 'sparseMatrix' class")
  } 
  if(!inherits(x, "dgCMatrix")){
    x <- as(x, "dgCMatrix") 
  }
  
  sdx <- rowSDs(x) 
  if(any(sdx==0)){
    zpr <- which(sdx==0)
    stop(paste0("Gene(s): "
                ,paste0(names(sdx)[zpr], collapse=", ")
                ,"  have standard deviation=0.  Please remove from them from counts matrix 'x'"
                )
         )
  } 
  
  return(x)
}

check_totalcounts <- function(counts_matrix, totalcounts){
  
  ### get totalcounts if not specified 
  if(is.null(totalcounts)){
    totalcounts <- Matrix::colSums(counts_matrix) 
  }
  stopifnot(length(totalcounts) == ncol(counts_matrix))
  if(!is.null(names(totalcounts)) && !is.null(colnames(counts_matrix))){
    stopifnot("ids of 'totalcounts' dont match ids in 'counts_matrix'" = 
                all(names(totalcounts) %in% colnames(counts_matrix))) 
    totalcounts <- totalcounts[colnames(counts_matrix)] 
  }
  if(is.null(names(totalcounts))) names(totalcounts) <- colnames(counts_matrix) 
  totalcounts <- pmax(totalcounts, 1)
  return(totalcounts) 
}

check_phi <- function(x, phi){
  stopifnot(length(phi) >=1)
  stopifnot(is.numeric(phi))
  if(length(phi) > 1){
    stopifnot(length(phi)==nrow(x))
    if(!is.null(names(phi)) && !is.null(rownames(x))){
      stopifnot("names of 'phi' dont match rownames in provided counts matrix 'x'" = 
                  all(names(phi) %in% rownames(x))) 
      phi <- phi[rownames(x)] 
    }
    if(is.null(names(phi))) names(phi) <- rownames(x) 
    stopifnot(all(phi > 0))
  }
  return(phi) 
}


check_grate <- function(x, grate, totalcounts){
  
  if(is.null(grate)){
    ### gene frequency relative to totalcounts, so subset counts matrices don't bias grate upward
    grate <- Matrix::rowSums(x) / sum(totalcounts) 
  }
  if(length(grate) != nrow(x)){
    stop(paste0("Mismatch between number of genes specified in gene frequency 'grate'\n"
                ,"and the number of rows in counts matrix 'x'.\n"
                ,"Was this intentional? If intending to use a subset of genes, \n"
                ,"please explicitly subset *both* the counts matrix 'x' and the gene frequencies 'grate'\n"
                ,"to matching dimensions to avoid silent errors or confusion.\n"
                ,"i.e., `x=counts[mygenes,], grate = genefreq[mygenes]"))
  }
  if(!is.null(names(grate)) && !is.null(rownames(x))){
    stopifnot("names of 'grate' dont match rownames in provided counts matrix 'x'" = 
                all(names(grate) %in% rownames(x))) 
    grate <- grate[rownames(x)] 
  }
  if(is.null(names(grate))) names(grate) <- rownames(x) 
  stopifnot("gene frequency must be >0 for all genes, otherwise gene should be removed from counts matrix or grate manually re-specified i.e. all(grate > 0)" = all(grate > 0))
  return(grate) 
}


check_grate_batch <- function(x, grate, batch_mat, totalcounts){
  
  if(is.null(grate)){
    ### expression rates by batch (\hat{p}), relative to per-batch totalcounts so subset counts matrices don't bias grate upward
    batch_totalcounts <- as.vector(totalcounts %*% batch_mat)
    grate <- x %*% batch_mat %*% Matrix::Diagonal(x = 1/batch_totalcounts, names = colnames(batch_mat))
  }
  stopifnot("For batch PCA, 'grate' should be a genes x batches matrix of frequencies" = inherits(grate, "Matrix"))
  stopifnot(ncol(grate)==ncol(batch_mat))
  stopifnot(colnames(grate)==colnames(batch_mat))
  if(nrow(grate) != nrow(x)){
    stop(paste0("Mismatch between number of genes specified in gene frequency 'grate'\n"
                ,"and the number of rows in counts matrix 'x'.\n"
                ,"Was this intentional? If intending to use a subset of genes, \n"
                ,"please explicitly subset *both* the counts matrix 'x' and the gene frequencies 'grate'\n"
                ,"to matching dimensions to avoid silent errors or confusion.\n"
                ,"i.e., `x=counts[mygenes,], grate = genefreq[mygenes]"))
  }
  if(!is.null(rownames(grate)) && !is.null(rownames(x))){
    stopifnot("rownames of 'grate' dont match rownames in provided counts matrix 'x'" = 
                all(rownames(grate) %in% rownames(x))) 
    grate <- grate[rownames(x),] 
  }
  if(is.null(rownames(grate))) rownames(grate) <- rownames(x) 
  stopifnot("gene frequency must be >0 for all genes in all batches, otherwise either gene should be removed or batch re-defined i.e. all(grate@x > 0)" = all(grate@x > 0))
  return(grate) 
}

