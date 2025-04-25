

#' Compute fold-change, proportion of cells >0, average expression metrics for each gene in each cluster 
#' (comparing with all other clusters)
#' 
#' @description
#' This function computes metrics by gene by cluster which can be passed as input to the `marker_heatmap()` function.
#' 
#' 
#' @param counts genes x cells matrix of counts.  this only needs to be provided if `normed` is not specified.  
#' When `normed` is not specified, `counts` will be 'totalcount-normalized' before 
#' foldchange metrics are calculated.
#' @param normed genes x cells matrix of data to be used for calculating fold change statistics by clusters.
#' If this is not provided, counts will be 'totalcount-normalized'.
#' If you desire to use raw counts for computing fold change, use like `normed=my_raw_counts`.
#' @param metadata metadata including (at minimum) 'cellid_column' and a 'cluster_column'
#' @param cluster_column column in metadata corresponding to cell type 
#' @param totalcounts optional user-specified vector of totalcounts used for normalizing the counts matrix.  
#' (useful if a subsetted counts matrix is passed)
#' 
#' @export
clusterwise_foldchange_metrics <- function(counts=NULL, normed = NULL, totalcounts = NULL, metadata, cluster_column){
  if(missing(normed)){
    normed <- Matrix::t(totalcount_norm(Matrix::t(counts), totalcounts))
  }
 
  metainfo <- data.table::copy(data.table::data.table(metadata))
  
  pb <- txtProgressBar(min =0, max = length(unique(metainfo[[cluster_column]])), style = 3)
  idx <- 0
  
  outl <- list()
  for(ii in unique(metainfo[[cluster_column]])){
    cells_ii <- metainfo[metainfo[[cluster_column]]==ii,cell_ID]
    cells_iiprime <- metainfo[metainfo[[cluster_column]]!=ii,cell_ID]
    cluster_expr_ii  <- Matrix::rowMeans(normed[,cells_ii,drop=FALSE])
    cluster_expr_iiprime <- Matrix::rowMeans(normed[,cells_iiprime,drop=FALSE]) 
    cluster_prop_ii <- Matrix::rowMeans(normed[,cells_ii,drop=FALSE] > 0) 
    cluster_prop_iiprime <- Matrix::rowMeans(normed[,cells_iiprime,drop=FALSE] > 0) 
    fctbl <- data.table::data.table(cluster=ii
                        ,cluster_expr = cluster_expr_ii
                        ,clusterprime_expr = cluster_expr_iiprime
                        ,gene = names(cluster_expr_ii)
                        ,cluster_prop = cluster_prop_ii
                        ,clusterprime_prop = cluster_prop_iiprime
                        ,ncells = length(cells_ii)
    )[,typ:=cluster_column]
    fctbl[,fold_change:=cluster_expr / clusterprime_expr]
    fctbl[,fold_change_prop:=cluster_prop / clusterprime_prop]
    outl[[paste0(ii)]] <- copy(fctbl) 
    idx <- idx + 1
    setTxtProgressBar(pb, idx)
  } 
  return(fc = data.table::rbindlist(outl))
}




