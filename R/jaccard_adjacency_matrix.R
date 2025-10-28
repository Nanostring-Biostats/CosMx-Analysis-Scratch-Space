

#' Create a pair of adjacency matrix graphs via jaccard similarity index.
#' 
#' @description
#' Analogous to Seurat::FindNeighbors, but computes across chunks of cells, and avoids CHOLMOD memory errors if numbers of cells is very large.
#' 
#' 
#' @param cell_embeddings cells x embeddings matrix (i.e., pca embeddings)  used to compute nearest neighbors.
#' @param k_nearest_neighbors number of nearest neighbors
#' @param annoy_metric  distance metric passed to `uwot:::annoy_nn`
#' @param npcs number of embedding dimensions to use from the cell_embeddings matrix (i.e., the number of PCs)
#' @param prune.SNN If the number of shared nearest neighbors is less than this value, trim the adjacency graph value to 0.
#' @param verbose Show the progress output from `annoy_nn`
#' 
#' @export
jaccard_adjacency_matrix <- function(cell_embeddings = NULL, k_nearest_neighbors = 20, annoy_metric = 'euclidean'
                                     ,npcs = NULL, prune.SNN = 1/15, verbose = TRUE){
  annoynn <- getFromNamespace("annoy_nn", "uwot")
  if(is.null(npcs)) npcs <- ncol(cell_embeddings)
  nnidx <-  
    annoynn(X = cell_embeddings[,1:npcs]
            ,metric = annoy_metric
            ,k = k_nearest_neighbors
            ,verbose = verbose
            ,search_k = -1
    )
  
  graphs <- list(nn = NULL, snn = NULL)
  graphs[["nn"]] <-  
    Matrix::sparseMatrix(i = rep(1:nrow(nnidx$idx), times = k_nearest_neighbors)
                         ,j = as.vector(nnidx$idx)
                         ,x = 1)
  rm(nnidx)
  gc()
  ncells <- ncol(graphs[["nn"]])
  
  ## compute the shared neighbors in smaller chunks, otherwise memory CHOLMOD errors if too many cells
  ### then get jaccard distance and prune
  chunkdt <- data.table::data.table(strt=unique(c(seq(1, ncells, 10**5), ncells)))
  chunkdt[,stp:=data.table::shift(strt, 1, type="lead")]
  chunkdt <- chunkdt[-.N]
  if(nrow(chunkdt) > 1){
    chunkdt[1:(.N-1),stp:=stp-1]
  }
  snnl <- vector(mode='list',length=nrow(chunkdt))
  
  message("Computing SNN") 
  pb <- txtProgressBar(min = 0, max = nrow(chunkdt), initial = 0,style=3) 
  for(chunk in 1:nrow(chunkdt)){
    strt <- chunkdt[chunk][["strt"]]
    stp <- chunkdt[chunk][["stp"]]
    snnl[[chunk]] <- graphs[["nn"]] %*% Matrix::t(graphs[["nn"]][strt:stp,])
    
    ### Jaccard calculation, intersection / (union)
    snnl[[chunk]]@x <- snnl[[chunk]]@x / (2*k_nearest_neighbors - snnl[[chunk]]@x)
    
    ### prune and drop explicit 0's  
    snnl[[chunk]]@x[snnl[[chunk]]@x < prune.SNN] <- 
      rep(0, sum(snnl[[chunk]]@x < prune.SNN))
    snnl[[chunk]] <-  Matrix::drop0(snnl[[chunk]])
    setTxtProgressBar(pb, chunk)
  }
  graphs[["snn"]] <- do.call(cbind, snnl) 
  rm(snnl)
  gc()
  dimnames(graphs[["nn"]]) <- dimnames(graphs[["snn"]]) <- list(c(rownames(cell_embeddings))
                                                                ,c(rownames(cell_embeddings))
  )
  return(graphs)
}


