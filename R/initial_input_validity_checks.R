

check_x_is_dgcmatrix <- function(x){
  if(!inherits(x, "sparseMatrix")){
    stop("'x' should be of 'sparseMatrix' class")
  } 
  if(!inherits(x, "dgCMatrix")){
    x <- as(x, "dgCMatrix") 
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
  return(totalcounts) 
}

check_grate <- function(x, grate){
  
  if(is.null(grate)){
    grate <- Matrix::rowSums(x) 
    grate <- grate / sum(grate) 
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
  return(grate) 
}


check_grate_batch <- function(x, grate, batch_mat){
  
  if(is.null(grate)){
    grate <- Matrix::rowSums(x) 
    grate <- grate / sum(grate) 
  }
  if(length(grate) != nrow(counts_matrix)){
    stop(paste0("Mismatch between number of genes specified in gene frequency 'grate'\n"
                ,"and the number of rows in counts matrix 'x'.\n"
                ,"Was this intentional? If intending to use a subset of genes, \n"
                ,"please explicitly subset *both* the counts matrix 'x' and the gene frequencies 'grate'\n"
                ,"to matching dimensions to avoid silent errors or confusion.\n"
                ,"i.e., `x=counts[mygenes,], grate = genefreq[mygenes]"))
  }
  if(!is.null(names(grate)) && !is.null(rownames(counts_matrix))){
    stopifnot("names of 'grate' dont match rownames in provided counts matrix 'x'" = 
                all(names(grate) %in% rownames(counts_matrix))) 
    grate <- grate[rownames(counts_matrix)] 
  }
  if(is.null(names(grate))) names(grate) <- rownames(counts_matrix) 
  return(grate) 
}