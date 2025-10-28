

#' Compute fold-change, proportion of positive cells, average expression metrics for each protein in each cluster 
#' (comparing with all other clusters)
#' 
#' @description
#' This function computes metrics by protein and by cluster which can be passed as input to the `marker_heatmap()` function.
#' 
#' 
#' @param raw proteins x cells expression matrix.  this only needs to be provided if `normed` is not specified.  
#' When `normed` is not specified, `raw` will be 'totalcount-normalized' before 
#' foldchange metrics are calculated.
#' @param normed proteins x cells matrix of data to be used for calculating fold change statistics by clusters.
#' If this is not provided, `raw` will be 'totalexpr-normalized'.
#' If you desire to use raw expression for computing fold change, use like `normed=my_raw_expr`.
#' @param metadata metadata including (at minimum) 'cellid_column' and a 'cluster_column'
#' @param cluster_column column in metadata corresponding to cell type 
#' @param cellid_column column in metadata corresponding to cell id. Should also correspond to the column names of the `raw` and `normed` expression matrices.
#' @param totalexpr optional user-specified cell-length vector of 'totalexpression' used for normalizing the raw matrix.  
#' (useful if a subsetted raw matrix is passed)
#' @param propd an expression matrix used for computing the proportion of 'positive cells' for a protein.  
#' For example, given a pearson-residual normalized expression matrix (`pearson_norm_mat`), we could pass 
#' `propd = (pearson_norm_mat > 0)`, or to use a stricter threshold,  `propd = (pearson_norm_mat > 1)`
#' 
#' @export
clusterwise_foldchange_metrics_protein <- 
  function(raw=NULL, normed = NULL, totalexpr = NULL, metadata, cluster_column, propd = NULL){
    
  stopifnot(cluster_column %in% colnames(metadata)) 
  metainfo <- data.table::copy(data.table::data.table(metadata))
  if(!(cellid_column) %in% colnames(metainfo)){
    stopifnot(!is.null(rownames(metadata)))
    cellid_column <- "cell_ID"
    metainfo[[cellid_column]] <- rownames(metadata)
  }
  if(cellid_column!="cell_ID" & ("cell_ID" %in% names(metainfo))) metainfo[["cell_ID"]] <- NULL
  data.table::setnames(metainfo, old=cellid_column, new="cell_ID")
  rm(metadata); gc()
  
  stopifnot("provided 'cellid_column' are not all unique in metadata " = 
              length(unique(metainfo[["cell_ID"]])) == nrow(metainfo))
  
  if(!is.null(raw)){
    ## checks on the raw matrix 
    if(!is.null(totalexpr)){
      stopifnot(length(totalexpr) == ncol(raw))
      if(!is.null(names(totalexpr)) && !is.null(colnames(raw))){
        stopifnot("ids of 'totalexpr' dont match ids in 'raw'" = 
                    all(names(totalexpr) %in% colnames(raw))) 
        totalexpr <- totalexpr[colnames(raw)] 
      }
      if(is.null(names(totalexpr))){
        warning("No names provided for 'totalexpr'.\nAssuming that 'totalexpr' vector matches the order of provided 'raw' matrix.\nPlease consider providing a cellid-named vector of totalexpr.")
      }
    }
    stopifnot(all(metainfo[["cell_ID"]] %in% colnames(raw)))
    if(nrow(metainfo)!=ncol(raw)){
      warning("Number of rows (cells) in 'metadata' does not match the number of columns (cells) in the 'raw' matrix!\nUsing only the cells in the metadata to compute the foldchange table.")
    }
  }  else {
    if(missing(normed)){
        stop("either 'normed' or 'raw' argument must be provided")
    }
  }
  
  if(missing(normed)){
    normed <- Matrix::t(totalcount_norm(Matrix::t(raw), totalexpr))
  }
 
  ## check on the normed matrix 
  stopifnot(all(metainfo[["cell_ID"]] %in% colnames(normed)))
  if(nrow(metainfo)!=ncol(normed)){
    warning("Number of rows (cells) in 'metadata' does not match the number of columns (cells) in the 'normed' matrix!\nUsing only the cells in the metadata to compute the foldchange table.")
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
    if(!is.null(propd)){
      cluster_prop_ii <- Matrix::rowMeans(propd[,cells_ii,drop=FALSE] > 0) 
      cluster_prop_iiprime <- Matrix::rowMeans(propd[,cells_iiprime,drop=FALSE] > 0) 
    } else {
      cluster_prop_ii <- Matrix::rowMeans(normed[,cells_ii,drop=FALSE] > 0) 
      cluster_prop_iiprime <- Matrix::rowMeans(normed[,cells_iiprime,drop=FALSE] > 0) 
    }
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