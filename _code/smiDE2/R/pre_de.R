
#' Calculate cell adjacencies and neighboring cell expression of genes before DE analysis.
#'
#' @param adjacencies_only avoid pre-computing expression matrices, only pre-compute the adjacencies between cells. 
#' @param metadata cell-level metadata corresponding to the cells in counts, with one row per cell.
#' Should contain the columns of cell_type_metadata_colname, 
#' cellid_colname, split_neighbors_by_colname, sdimx_colname, and sdimy_colname.
#' @param cell_type_metadata_colname column (commonly cell type column) corresponding to categories by which cell expression is 
#' aggregated by to calculate expression among neighboring cells by cell_type_metadata_colname type.
#' @param cellid_colname column name corresponding to cell id.
#' @param sdimx_colname column name corresponding to 'x' axis spatial coordinate of cell.
#' @param sdimy_colname column name corresponding to 'y' axis spatial coordinate of cell.
#' @param split_neighbors_by_colname (optional) If specified, identification of neighboring cells is split by this 
#' column in the meta data.  For example, to avoid identifying cells as neighbors if they were to have close x/y coordinates, but were on different tissue.
#' @param mm_radius maximum euclidean distance in units of sdimx, sdimy by which two cells will be identified as neighbors. 
#' @param verbose print some extra messages during calculation.
#' @param counts matrix of counts, only required if adjacencies_only = FALSE
#' @param normalized_data optional, matrix of normalized data. 
#' If not provided, totalCounts normalization is applied to counts. Only relevant if adjacencies_only = FALSE
#' @param ref_celltype Reference category (commonly a cell type), used to calculate neighbor expression with respect to.
#' Can be vector of one or more levels of cell_type_metadata_colname. 
#' Can use keyword "all" for all levels of cell_type_metadata_colname.
#' Only relevant if adjacencies_only = FALSE
#' @param nb_celltypes_to_individually_calc vector of cell types for which to individually calculate neighbor expression.
#' default is NULL, which will make a separate list in neighbor_expr_by_ct output with matrix of expression by cell type.
#' Only relevant if adjacencies_only = FALSE
#' @param aggregation "sum" (default) calculates a total amount of target expressed in neighboring cells.  Anything other values passed will calculate an average.
#' Only relevant if adjacencies_only = FALSE
#' @param weight_colname (optional).  If NULL, neighbor expression is unweighted.  
#' If "weight", the column corresponding to "weight" (1/distance) by default is used to weight calculated expression of neighbor cells.
#' Only relevant if adjacencies_only = FALSE
#' 
#' @return list object with data.table of cell_adjacencies (one row per pair of adjacent cells),
#'  a list (nblist) with sparse matrices corresponding to neighbor cell expression passed to downstream DE in smi_de or run_de functions.
#'  
#' @examples 
#' library(Seurat)
#' library(data.table); setDTthreads(1)
#' datadir<-system.file("extdata", package="smiDE")
#' sem <- readRDS(paste0(datadir, "/small_nsclc.rds"))
#' metainfo <- data.table::data.table(sem@meta.data)
#'  
#' pre_de_obj <- 
#'    pre_de(metadata = metainfo
#'           ,adjacencies_only = TRUE
#'           ,cell_type_metadata_colname = "cell_type"
#'           ,split_neighbors_by_colname = "tissue"
#'           ,mm_radius = 0.05
#'           ,sdimx_colname = "sdimx"
#'           ,sdimy_colname = "sdimy"
#'           ,verbose=TRUE
#'    ) 
#'
#' # data.table of cell-cell adjacencies
#' pre_de_obj$cell_adjacency_dt[]
#' 
#' ### Note:  
#' ### The examples below create genes x cells neighbor-expression matrices.
#' ### This method for pre-calculating the large neighbor expression matrices
#' ### is high-memory and currently deprecated, in favor of doing these calculations
#' ### 'on-the-fly' when running smi_de.
#' ### However, the examples below will still run, and are useful for visual inspection
#' ### of what's being calculated 'under-the-hood' when running DE models and using
#' ### covariate-adjustment to control for neighbor expression.
#' 
#' totalcounts <- sem@meta.data[["totalcounts"]]
#' sem <- Seurat::SetAssayData(sem
#'                            ,"data"
#'                            ,sem[["RNA"]]@counts %*% Matrix::Diagonal(x=mean(totalcounts)/totalcounts, names=colnames(sem))
#'                            )
#' 
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
#' names(pre_de_obj)
#' 
#' # Example for measuring neighbor expression with 'all' levels of
#' # ref_celltype = 'all'
#' pre_de_obj <- 
#' pre_de(counts = sem[["RNA"]]@counts
#'        ,normalized_data = sem[["RNA"]]@data
#'        ,metadata = metainfo
#'        ,cell_type_metadata_colname = "cell_type"
#'        ,split_neighbors_by_colname = "tissue"
#'        ,mm_radius = 0.05
#'        ,ref_celltype = c("all")
#'        ,sdimx_colname = "sdimx"
#'        ,sdimy_colname = "sdimy"
#'        ,weight_colname = "weight"
#'        ,aggregation = "sum"
#'        ,verbose=TRUE
#' )
#' names(pre_de_obj)
#' dim(pre_de_obj$nblist$neighbor_expr_byct$otherct)  
#' metainfo[, .N]
#'
#'
#' # Example for measuring neighbor expression with a vector of multiple cell types
#' pre_de_obj <- 
#' pre_de(counts = sem[["RNA"]]@counts
#'        ,normalized_data = sem[["RNA"]]@data
#'        ,metadata = metainfo
#'        ,cell_type_metadata_colname = "cell_type"
#'        ,split_neighbors_by_colname = "tissue"
#'        ,mm_radius = 0.05
#'        ,ref_celltype = c("macrophage", "fibroblast")
#'        ,sdimx_colname = "sdimx"
#'        ,sdimy_colname = "sdimy"
#'        ,weight_colname = "weight"
#'        ,aggregation = "sum"
#'        ,verbose=TRUE
#' )
#' names(pre_de_obj)
#' dim(pre_de_obj$nblist$neighbor_expr_byct$otherct)  
#' metainfo[cell_type %in% pre_de_obj$nblist$ref_celltype, .N]
#'  
#' @export  
pre_de <- function(adjacencies_only=TRUE
                   ,metadata
                   ,ref_celltype
                   ,cell_type_metadata_colname
                   ,weight_colname = NULL
                   ,cellid_colname = "cell_ID"
                   ,sdimx_colname = "sdimx"
                   ,sdimy_colname = "sdimy"
                   ,split_neighbors_by_colname = "tissue" 
                   ,mm_radius = 0.05
                   ,counts = NULL
                   ,normalized_data = NULL
                   ,verbose=TRUE
                   ,aggregation = "sum"
                   ,nb_celltypes_to_individually_calc = NULL
                   ){

  metainfo <- data.table::copy(as.data.table(metadata))

  if(!adjacencies_only){
    if(missing(counts)) stop("counts argument must be provided if adjacencies_only=FALSE")
  }
  if(missing(cell_type_metadata_colname)) stop("cell_type_metadata_colname missing with no default.")
  stopifnot(cell_type_metadata_colname %in% colnames(metainfo))
  stopifnot(!is.null(rownames(metadata)) & (cellid_colname %in% names(metainfo)))
  if(!(cellid_colname) %in% names(metainfo)){
    cellid_colname <- "cell_ID"
    metainfo[[cellid_colname]] <- rownames(metadata)
  }
   
  stopifnot(sdimx_colname %in% names(metainfo)) 
  stopifnot(sdimy_colname %in% names(metainfo)) 
  
  if(!adjacencies_only){
    stopifnot(cell_type_metadata_colname %in% names(metainfo)) 
  } 
 
 
  if(!adjacencies_only){
    if(length(setdiff(colnames(counts), metainfo[[cellid_colname]])) > 0 ||
       length(setdiff(metainfo[[cellid_colname]], colnames(counts))) > 0 
       ){
      warning("Not identical set of cells between assay matrix and meta.data")
      mcells <- nrow(metainfo)
      acells <- ncol(counts)
      comm <- length(intersect(colnames(counts), metainfo[[cellid_colname]]))
      message(paste0(mcells, " cells in metadata."))
      message(paste0(acells, " cells in assay."))
      message(paste0(comm, " cells common between metadata and assay."))
    }
    if("all" %in% ref_celltype) ref_celltype <- unique(metainfo[[cell_type_metadata_colname]])
    stopifnot(all(ref_celltype %in% unique(metainfo[[cell_type_metadata_colname]])))
    stopifnot(all(setdiff(nb_celltypes_to_individually_calc, c("otherct", "allct")) %in% unique(metainfo[[cell_type_metadata_colname]])))
    
    if(is.null(normalized_data)){
      colsumms <- Matrix::colSums(counts)
      norm_factors <- mean(colsumms)/colsumms
      norm_factors[colsumms==0] <- 1
      normalized_data <- counts %*% Matrix::Diagonal(x=norm_factors, names=colnames(counts))
      #normalized_data <- Matrix::t(Matrix::t(counts)/(norm_factors))
      ### apply normalization factor by column (cell).
    }
    stopifnot(all.equal(dim(counts), dim(normalized_data)))
  } 
   
  if(is.null(split_neighbors_by_colname)){
    split_neighbors_by_colname <- "all_data"
    metainfo[[split_neighbors_by_colname]] <- "all_cells" 
  }
   
  # Identify cell-cell spatial neighbors 
  if(verbose) message(paste0(Sys.time()
                             ,", identifying cell-cell spatial neighbors within "
                             ,mm_radius, " radius."))
  
  tiss <- split(metainfo[,c(cellid_colname, sdimx_colname, sdimy_colname,split_neighbors_by_colname),with=FALSE]
                ,by=split_neighbors_by_colname)
  cell_adjacency_dt <-
    rbindlist(
      lapply(1:length(tiss), function(x){
        nbr <- fast_make_all_neighbors(tiss[[x]]
                                       ,cellid_col = cellid_colname
                                       ,sdimx_col = sdimx_colname
                                       ,sdimy_col = sdimy_colname
                                       ,radius = mm_radius
                                       )
        if(verbose) message(paste0("neighbors calculated for "
                                 ,split_neighbors_by_colname
                                 , ": "
                                 ,tiss[[x]][[split_neighbors_by_colname]][1]
                                 ))
        return(nbr)
      })
    ) 
  rm(tiss); gc()
  
  if(length(cell_type_metadata_colname) > 0){
    cell_adjacency_dt <- 
    merge(cell_adjacency_dt, metainfo[,c(cellid_colname, cell_type_metadata_colname),with=FALSE]
          ,by.x=c("from"), by.y=c(cellid_colname)
          ,sort=FALSE)
    cell_adjacency_dt <- 
    merge(cell_adjacency_dt, metainfo[,c(cellid_colname, cell_type_metadata_colname),with=FALSE]
          ,by.x=c("to"), by.y=c(cellid_colname)
          ,sort=FALSE
          ,suffixes=c("_from", "_to"))
  } 
  nblist <- NULL
  if(!adjacencies_only){
    
    # Create neighborhood list objects from normalized data. 
    if(verbose) message(paste0(Sys.time()
                               ,". creating neighborhood list object"
                               #, ref_celltype
                               #, " ."
                               ))
   
    ict <- 1L 
    if(verbose) message(paste0("Measuring neighbor expression for ", ref_celltype[1]))
    if(verbose) message(paste0("(cluster ",1, " of ", length(ref_celltype), ")"))
    nblist <- 
      measure_neighbor_expression_by_celltype(assay_matrix = normalized_data
                                              ,metadata = metainfo
                                              ,cell_adjacencies = cell_adjacency_dt
                                              ,ref_celltype[1]
                                              ,cell_type_metadata_colname
                                              ,aggregation
                                              ,cellid_colname
                                              ,weight_colname
                                              #,Wmat
                                              #,ct_matrix
                                              #,adjacency_counts_by_ct
                                              ,nb_celltypes_to_individually_calc = nb_celltypes_to_individually_calc
      )
    
    if(length(ref_celltype) > 1){
      for(ct in ref_celltype[2:length(ref_celltype)]){
        ict <- ict + 1L
        if(verbose) message(paste0(Sys.time()))
        if(verbose) message(paste0("Measuring neighbor expression for ", ct))
        if(verbose) message(paste0("(cluster ",ict, " of ", length(ref_celltype), ")"))
        #nblistl[[paste0(ct)]] <- 
        nblist_tmp <- 
          measure_neighbor_expression_by_celltype(assay_matrix = normalized_data
                                                  ,metadata = metainfo
                                                  ,cell_adjacencies = cell_adjacency_dt
                                                  ,ct
                                                  ,cell_type_metadata_colname
                                                  ,aggregation
                                                  ,cellid_colname
                                                  ,weight_colname
                                                  ,Wmat = nblist$adjacency_mat 
                                                  ,ct_matrix = nblist$ct_matrix
                                                  ,adjacency_counts_by_ct = nblist$adjacency_counts_by_ct 
                                                  ,nb_celltypes_to_individually_calc = nb_celltypes_to_individually_calc
          )
        
       #for(nm in names(nblist[["neighbor_expr_byct"]])){
       for(nm in names(nblist[["neighbor_expr_byct"]])){
         nblist[["neighbor_expr_byct"]][[nm]] <- 
           Matrix::cbind2(nblist[["neighbor_expr_byct"]][[nm]]
                          ,nblist_tmp[["neighbor_expr_byct"]][[nm]]
         )
       }
       rm(nblist_tmp); gc()
      }
    }
    
    nblist$ref_celltype <- ref_celltype 
    
  } 
  return_obj <- list(nblist = nblist
                    ,cell_adjacency_dt = cell_adjacency_dt
                    )
  class(return_obj) <- append(class(return_obj), "prede")
  return(return_obj)
}
  