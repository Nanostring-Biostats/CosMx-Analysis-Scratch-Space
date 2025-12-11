#' Clusters by gene overlap ratio metric.
#'
#' For each cluster and gene, returns 
#' * the average expression of `assay_matrix` corresponding to cells in a cluster
#' * the average expression of `assay_matrix` in cells which are both neighbors to cells in the cluster, and not belonging to the same cluster.
#'
#' Genes which are lowly expressed in a cluster, but highly expressed in **neighboring cells** of that cluster may be more prone to expression primarily from overlapping cells or segmentation errors. 
#' Overlap-prone genes may be prioritized for exclusion from downstream cell type specific analyses , such as cell type specific differential expression.
#' 
#' @param assay_matrix counts or normalized assay matrix.  
#'                     Raw counts recommended for count distributions s.a. negative binomial (i.e., family=nbinom2). 
#'                     Normalized is recommended for gaussian models.
#' @param metadata data.table or data.frame of meta.data associated with assay.  
#'                 Should include columns corresponding to "cluster_col" and "grouping_col" argument.  If `adjacency_matrix` is
#' @param adjacency_matrix cell x cell matrix indicating which cells are neighbors of which other cells.  
#' This can be passed as the nblist$cell_adjacency_mat object from pre_de object obtained from running smiDE::pre_de. 
#' If missing or NULL, this will be created from `sdimx_col`, `sdimy_col`, `cellid_col`
#' ,and optional `split_neighbors_by_colnames` columns in the metadata.
#' @param cluster_col column indicating clusters by which to calculate overlap metric. For example, column indicating cell type.
#' @param grouping_col optional vector of column names in meta data by which to separately calculate overlap metric. For example, could be the primary variable of interest in a Differential Expression (DE) analysis. 
#' For example, if running a DE analysis for 'niche', `grouping_col='niche'` will calculate metric separately for each 'niche', and can be used to define a filter for excluding genes likely to be 
#' expressed in overlapping cells for certain cell types in any level of 'niche'. 
#' @param cellid_col name of column indicating the cell ID in the metadata, i.e., "cell_ID", corresponding to column names in provided `assay_matrix`.
#' @param sdimx_col name of column with x-coordinate for cell, i.e., "sdimx". Only needed if adjacency_matrix is not provided.
#' @param sdimy_col name of column with x-coordinate for cell, i.e., "sdimy". Only needed if adjacency_matrix is not provided.
#' @param radius radius within which to identify all other neighboring cells, i.e, 0.05.  Only needed if adjacency_matrix is not provided.
#' @param verbose print message to screen describing status of calculation.
#' 
#' @examples
#' library(Seurat)
#' library(data.table); setDTthreads(1)
#' datadir<-system.file("extdata", package="smiDE")
#' sem <- readRDS(paste0(datadir, "/small_nsclc.rds"))
#'
#' totalcounts <- sem@meta.data[["totalcounts"]]
#' sem <- Seurat::SetAssayData(sem
#'                            ,"data"
#'                            ,sem[["RNA"]]@counts %*% Matrix::Diagonal(x=mean(totalcounts)/totalcounts, names=colnames(sem))
#'                            )
#'
#' 
#' metainfo <- data.table::data.table(sem@meta.data)
#' metric_by_niche_and_tissue <- 
#' overlap_ratio_metric(sem[["RNA"]]@data
#'                      ,metadata = metainfo
#'                      ,cluster_col = "cell_type"
#'                      ,grouping_col = c("niche", "tissue") ### calculate within each `niche` and each `tissue`
#'                      ,split_neighbors_by_colname = c("Run_Tissue_name")
#'                      ,radius = 0.05
#'                      ,sdimx_col = "sdimx"
#'                      ,sdimy_col = "sdimy"
#'                      ,cellid_col = "cell_ID"
#'                      )
#'                            
#' ### top marker in 'macrophages' within each niche and tissue category
#' metric_by_niche_and_tissue[cell_type=="macrophage"][order(ratio),head(.SD,1),by=.(niche,tissue)]
#' 
#' ### top expressor in 'other cell types' within each niche and tissue category
#' metric_by_niche_and_tissue[cell_type=="macrophage"][order(-ratio),head(.SD,1),by=.(niche,tissue)]
#'
#' 
#' metric_overall_1 <- 
#' overlap_ratio_metric(sem[["RNA"]]@data
#'                      ,metadata = metainfo
#'                      ,cluster_col = "cell_type"
#'                      ,grouping_col = NULL
#'                      ,split_neighbors_by_colname = c("Run_Tissue_name")
#'                      ,radius = 0.05
#'                      ,sdimx_col = "sdimx"
#'                      ,sdimy_col = "sdimy"
#'                      ,cellid_col = "cell_ID"
#'                      )
#' metric_overall_1[cell_type=="macrophage"][order(ratio)][1:10] ## highly expressed in macrophages
#' metric_overall_1[cell_type=="macrophage"][order(-ratio)][1:10] ## highly expressed in macrophage neighbors
#' 
#' ### Using a adjacency matrix from metadata. (Haven't run smiDE::pre_de)
#' ### Using adjacency matrix from a pre-calculated pre_de object with adjacencies_only = FALSE
#' pre_de_obj <- 
#' pre_de(adjacencies_only = FALSE
#'        ,counts = sem[["RNA"]]@counts
#'        ,normalized_data = sem[["RNA"]]@data
#'        ,metadata = metainfo
#'        ,cell_type_metadata_colname = "cell_type"
#'        ,split_neighbors_by_colname = "tissue"
#'        ,mm_radius = 0.05
#'        ,ref_celltype = "fibroblast"
#'        ,sdimx_colname = "sdimx"
#'        ,sdimy_colname = "sdimy"
#'        ,weight_colname = "weight"
#'        ,aggregation = "sum"
#'        ,verbose=TRUE
#'  ) 
#' 
#' ### Using adjacency matrix from a pre-calculated pre_de object with adjacencies_only = FALSE
#' metric_overall_2 <- 
#' overlap_ratio_metric(sem[["RNA"]]@data
#'                      ,metadata = sem@meta.data
#'                      ,cluster_col = "cell_type"
#'                      ,grouping_col = NULL
#'                      ,adjacency_matrix = pre_de_obj$nblist$adjacency_mat
#'                      ,cellid_col = "cell_ID"
#'                      )
#' 
#' 
#' metric_overall_2[1:10]
#' 
#' metric_overall_2[1:5]; metric_overall_1[1:5] ## same
#' 
#' @export
#' 
#' 

overlap_ratio_metric <- function(assay_matrix
                                 ,metadata
                                 ,adjacency_matrix = NULL
                                 ,cluster_col = c("RNA_nbclust_clusters", "spatialClusteringAssignments")
                                 ,cellid_col = "cell_ID"
                                 ,grouping_col = NULL 
                                 ,sdimx_col = NULL
                                 ,sdimy_col = NULL
                                 ,radius = NULL
                                 ,split_neighbors_by_colname = NULL
                                 ,verbose = FALSE
                                 ){
 
  metainfo <- data.table::copy(as.data.table(metadata))
  stopifnot(all(c(cluster_col, cellid_col, grouping_col, split_neighbors_by_colname) %in% names(metainfo)))
   
  if(is.null(adjacency_matrix)){
    if(missing(sdimx_col)) stop("if adjacency matrix not provided, must provide x coordinate column 'sdimx_col'.")
    if(missing(sdimy_col)) stop("if adjacency matrix not provided, must provide y coordinate column 'sdimy_col'.")
    if(missing(radius)) stop("if adjacency matrix not provided, must provide maximum distance for identifying neighbors via 'radius' argument.")
    
    if(is.null(split_neighbors_by_colname)){
      split_neighbors_by_colname <- "all_data"
      metainfo[[split_neighbors_by_colname]] <- "all_cells" 
    }
    stopifnot(all(c(sdimx_col, sdimy_col, split_neighbors_by_colname) %in% colnames(metainfo)))
    
    tiss <- split(metainfo[,c(cellid_col, sdimx_col, sdimy_col,split_neighbors_by_colname),with=FALSE]
                  ,by=split_neighbors_by_colname)
    cell_adjacency_dt <-
      rbindlist(
        lapply(1:length(tiss), function(x){
          nbr <- fast_make_all_neighbors(tiss[[x]]
                                         ,cellid_col
                                         ,sdimx_col
                                         ,sdimy_col
                                         ,radius
          )
          if(verbose) message(paste0("neighbors calculated for "
                                     ,paste0(split_neighbors_by_colname, collapse=".")
                                     , ": "
                                     ,names(tiss)[x]
                                     )
          )
          return(nbr)
        })
      ) 
    rm(tiss); gc() 
    W <- smiDE:::make_W(metainfo
                   ,cell_adjacency_dt
                   ,cluster_col
                   ,aggregation = 'sum' 
                   ,cellid_colname = cellid_col
                   ,weight_colname = NULL
    )
  } else {
    W <- adjacency_matrix 
  }

  W@x[W@x > 0] <- 1 ## override non-zero or distance based weights in adjacency matrix and 
                    ## replace with 1's, used to take a simple average of expression thats interpretable and comparable 
                    ## between cluster and neighbors.
  
  if(!Matrix::isSymmetric(W)) warning("provided adjacency matrix is not symmetric") 
  stopifnot(all(colnames(assay_matrix) %in% rownames(W)))
  stopifnot(all(rownames(W) == colnames(W)))
  
  W <- W[colnames(assay_matrix),colnames(assay_matrix)]
   
  if(is.null(grouping_col)){
    grouping_col <- "all_data"
    metainfo[[grouping_col]] <- "all_cells"
  }
  grpsplit <- split(metainfo, by = grouping_col)
  metricl <- list() 
  for(jj in 1:length(grpsplit)){
    metricl_grp <- list()
    if(verbose) message(names(grpsplit)[jj])
    xx <- grpsplit[[jj]]
    grpcells <- xx[[cellid_col]]
    cellsplit <- split(xx, by = cluster_col)
    for(ii in 1:length(cellsplit)){
      if(verbose) message(names(cellsplit)[ii])
      cells_ii <- cellsplit[[ii]][[cellid_col]]
      other_cells <- setdiff(grpcells, cells_ii)
      
      ## column standardized in order to take the average
      rs <- Matrix::colSums(W[other_cells,,drop=FALSE])  ## denominator for number of cells of other cell types
      rs[rs==0] <- Inf
      M <- Matrix::sparseMatrix(i=1:nrow(W), j=1:ncol(W),x=1/rs 
                                ,dimnames=dimnames(W))
      
      Wcs <- W %*% M ## column standardized adjacency matrix.  (columns sum to 1.)
      
      ct_mean_expr <- Matrix::rowMeans(assay_matrix[,cells_ii,drop=FALSE]) ## average expression of group
      
      ## (genes x other cells) * (other cells x cells in group)
      other_ii <- assay_matrix[,other_cells,drop=FALSE] %*% 
        Wcs[other_cells,cells_ii,drop=FALSE] 
      
      avg_percell <- Matrix::rowMeans(other_ii) ## average expression of neighboring cells of group
      avg_percell <- data.table(target=names(avg_percell), avg_neighbor_othercluster = avg_percell)
      ct_mean_expr <- data.table(target=names(ct_mean_expr), avg_cluster = ct_mean_expr)
      metricl_grp[[ii]] <- merge(avg_percell
                           ,ct_mean_expr
                           ,by="target")
     
      metricl_grp[[ii]][[cluster_col]] <- names(cellsplit)[ii] ## add annotation for cluster column
      for(kk in grouping_col){ ## add annotation for grouping_col
        metricl_grp[[ii]][[kk]] <- cellsplit[[ii]][1][[kk]] 
      }
    } 
    metricl[[jj]] <-rbindlist(metricl_grp)
  } 
   
  pdmetric <- rbindlist(metricl)  
  pdmetric[,ratio:=avg_neighbor_othercluster/avg_cluster]
  return(pdmetric)
}

