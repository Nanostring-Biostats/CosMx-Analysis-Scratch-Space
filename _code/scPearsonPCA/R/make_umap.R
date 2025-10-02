

#' A UMAP helper function 
#' 
#' @description
#' Create a UMAP from a PCA object and also return the fuzzy nearest-neighbors graph used to 
#' build the UMAP.  This nearest neighbor graph can be used as the starting point for unsupervised clustering, 
#' potentially improving the agreement between UMAP and unsupervised clusters derived from the same PC data.
#'
#' @param pcaobj A pca list object containing a DimReduc object called 'reduction.data' which is passed to `uwot::umap`. Same format as that returned by `sparse_quasipoisson_pca_seurat`.  
#' @param min_dist  hyperparameter passed directly to `uwot::umap`
#' @param n_neighbors  hyperparameter passed directly to `uwot::umap`
#' @param metric  hyperparameter passed directly to `uwot::umap`
#' @param key  DimReduc key for umap object.
#' subset (of, say, highly variable) genes to be used in PCA. 
#'
#'
#' @export  
make_umap <- function(pcaobj, min_dist=0.01, n_neighbors=30, metric="cosine",key ="UMAP_" ){
  ump <- 
    uwot::umap(pcaobj$reduction.data@cell.embeddings
               ,n_neighbors = n_neighbors
               ,nn_method = "annoy"
               ,metric = metric
               ,min_dist = min_dist
               ,ret_extra = c("fgraph","nn")
               ,verbose = TRUE)
  
  umpgraph <- ump$fgraph
  dimnames(umpgraph) <- list(rownames(ump$nn[[1]]$idx), rownames(ump$nn[[1]]$idx))
  colnames(ump$embedding) <- paste0(key, c(1,2)) 
  ump <- Seurat::CreateDimReducObject(embeddings = ump$embedding, key = key)
  return(list(grph = umpgraph
              ,ump = ump))
}

