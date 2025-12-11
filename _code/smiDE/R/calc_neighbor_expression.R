
make_ct_matrix <- function(
                           metadata
                           ,cell_adjacencies
                           ,cell_type_metadata_colname
                           ,cellid_colname="cell_ID"
  ){
  metainfo <- data.table::copy(as.data.table(metadata))
  stopifnot(cell_type_metadata_colname %in% names(metainfo)) 
  
  setnames(metainfo
          ,old=c(cell_type_metadata_colname,cellid_colname)
          ,new=c("cell_type_metadata_colname__"
                 ,"cell_ID"))
  
  
  ct_dt <- metainfo[,.(cell_ID, cell_type_metadata_colname__)] 
  setkeyv(ct_dt, "cell_type_metadata_colname__")
  
  ct_dt[,jj:=.GRP,by=cell_type_metadata_colname__] 
  ct_dt <- ct_dt[match(metainfo$cell_ID,cell_ID)]
  ct_dt[,ii:=.GRP,by=cell_ID]
  
  ct_matrix <- Matrix::sparseMatrix(i=ct_dt$ii
                                    ,j=ct_dt$jj
                                    ,x=1L
                                    ,giveCsparse = TRUE
                                    ,dimnames = 
                                      list(c(ct_dt[order(ii),unique(cell_ID)])
                                           ,c(ct_dt[order(jj),unique(cell_type_metadata_colname__)])) 
                                    ) 
  return(ct_matrix)
}

make_W <- function(
                   metadata
                   ,cell_adjacencies
                   ,cell_type_metadata_colname
                   ,aggregation="sum"
                   ,cellid_colname="cell_ID"
                   ,weight_colname=NULL
                   ,standardize_matrix_by_row_or_col = c("col", "row")
  ){

  metainfo <- data.table::copy(data.table::as.data.table(metadata)) 
  if(!(cellid_colname) %in% names(metainfo)){
    cellid_colname <- "cell_ID"
    metainfo[[cellid_colname]] <- rownames(metadata)
  }
  if(cellid_colname!="cell_ID" & "cell_ID" %in% names(metainfo)) metainfo[["cell_ID"]] <- NULL
  setnames(metainfo, old=cellid_colname, new="cell_ID")
  
  standardize_matrix_by_row_or_col <- match.arg(standardize_matrix_by_row_or_col) 
  idx <- rbindlist(list(cell_adjacencies[,.(cell_ID=from)][,unique(.SD)]
                        ,cell_adjacencies[,.(cell_ID=to)][,unique(.SD)]
                       )
                   )[,unique(.SD)]
  
  ### Add in any missing cells with 0 adjacencies
  idx <- rbindlist(list(idx
                        ,data.table(cell_ID=setdiff(metainfo[,unique(cell_ID)]
                                                    ,idx[,cell_ID]
                       )
                  )
  ))
  idx <- idx[order(cell_ID)]
  idx[,ind:=1:.N]
   
   
  ### Join numeric column/row index onto cell_id from and to 
  cell_adjacencies <- merge(cell_adjacencies, idx[,.(from=cell_ID, ind1=ind)], by=c("from"))
  cell_adjacencies <- merge(cell_adjacencies, idx[,.(to=cell_ID, ind2=ind)], by=c("to"))
  
  if(is.null(weight_colname)){
    cell_adjacencies[,wght__:=1L]
  }  else {
    cell_adjacencies[,wght__:=cell_adjacencies[[weight_colname]]] 
  }
  
  if(!idx[,.N] %in% cell_adjacencies[,unique(c(ind1, ind2))]){
    fillind <- data.table(ind1=1,ind2=idx[,.N],wght__=0)
    cell_adjacencies <- rbindlist(list(cell_adjacencies
                                           ,fillind),use.names=TRUE,fill=TRUE)
  }
  ### Make sparse neighborhood adjacency matrix.  
  ### each row is an adjacency belonging to row "i" and column "j"
  ### , as well as row "j" and column "i" for symmetry.
  Wmat <- Matrix::sparseMatrix(i=c(cell_adjacencies[,ind1], cell_adjacencies[,ind2])
                               ,j=c(cell_adjacencies[,ind2],cell_adjacencies[,ind1])
                               ,x=c(cell_adjacencies[,wght__],cell_adjacencies[,wght__])
                               ,dimnames = list(idx[,cell_ID],idx[,cell_ID])
                               ,repr="T")
  ###  
  if(aggregation=="sum"){
    Mmat <- Matrix::Diagonal(idx[,.N])
    dimnames(Mmat) <- list(idx[,cell_ID],idx[,cell_ID])
  } else {
    rs <- Matrix::colSums(Wmat) 
    rs[rs==0] <- Inf
    Mmat <- Matrix::sparseMatrix(i=1:idx[,.N], j=1:idx[,.N],x=1/rs 
                                 ,dimnames=list(idx[,cell_ID],idx[,cell_ID])
                                 ,repr="T")
  }
  if(standardize_matrix_by_row_or_col=="row"){
    Wmat <- Mmat %*% Wmat
  } else {
    ### Column standardized W matrix.  
    Wmat <- Wmat %*% Mmat
  }
  
  return(Wmat)
}

#' Measure one gene's expression in neighboring cells, on-the-fly in DE
#'
measure_neighbor_expression_single <- 
  function(assay_matrix
           ,metadata
           ,cell_adjacencies
           ,ref_celltype
           ,cell_type_metadata_colname
           ,aggregation="sum"
           ,cellid_colname="cell_ID"
           ,weight_colname=NULL
           ,Wmat=NULL
#           ,ct_matrix = NULL
           ,adjacency_counts_by_ct = NULL
           ,nb_celltypes_to_individually_calc = NULL
  ){
    
#    metainfo <- data.table::copy(as.data.table(metadata))
#    
#    stopifnot(cell_type_metadata_colname %in% names(metainfo)) 
#    stopifnot(ref_celltype %in% unique(metainfo[[cell_type_metadata_colname]]))
#    
#    setnames(metainfo
#             , old=c(cell_type_metadata_colname,cellid_colname)
#             ,new=c("cell_type_metadata_colname__"
#                    ,"cell_ID"))
#
#    
#    if(missing(ct_matrix)){
#      ct_matrix <- make_ct_matrix(metadata
#                                  ,cell_adjacencies
#                                  ,cell_type_metadata_colname
#                                  ,cellid_colname
#      )
#    }   
    if(missing(Wmat)){
      Wmat <- make_W(metadata
                     ,cell_adjacencies
                     ,cell_type_metadata_colname
                     ,aggregation
                     ,cellid_colname
                     ,weight_colname
                     )
    } 
      
     
    # extract expression matrix to use, align columns with rows of adjacency matrix for matrix multiplication
    ex <- assay_matrix[,rownames(Wmat),drop=FALSE]
    rm(assay_matrix); gc()
    
    ct_ref <- metainfo[cell_type_metadata_colname__==ref_celltype,cell_ID]
   
    uniq_ct <- metainfo[,unique(cell_type_metadata_colname__)]
    if(is.null(nb_celltypes_to_individually_calc)){
      nb_celltypes_to_individually_calc <- uniq_ct
      nb_celltypes_to_individually_calc <- c("allct", "otherct", nb_celltypes_to_individually_calc) 
    } 
    ind_ct_to_calc <- setdiff(nb_celltypes_to_individually_calc, c("otherct", "allct"))
    stopifnot(all(ind_ct_to_calc %in% uniq_ct))

    ## Neighboring counts by individual cell types 
    neighbor_expr_byct <- vector(mode="list", length = length(nb_celltypes_to_individually_calc))
    names(neighbor_expr_byct) <- gsub("\\-|\\ ","_",nb_celltypes_to_individually_calc)
    
    if("allct" %in% nb_celltypes_to_individually_calc){
      ## All neighboring counts
      neighbor_expr_byct[["allct"]] <- 
        ex %*%
        Wmat[,ct_ref,drop=FALSE]
    }
    
    if("otherct" %in% nb_celltypes_to_individually_calc){
      ## All neighboring counts except those of the same cell type
      neighbor_expr_byct[["otherct"]] <- 
        ex[,setdiff(colnames(ex),ct_ref),drop=FALSE] %*%
        Wmat[setdiff(colnames(ex),ct_ref),ct_ref,drop=FALSE]
      
    } 
    
    if(length(ind_ct_to_calc) > 0){
      for(xx in ind_ct_to_calc){
        ctidx <- metainfo[cell_type_metadata_colname__==xx,cell_ID] 
        neighbor_expr_byct[[gsub("\\-|\\ ","_",xx)]] <-   
          ex[,ctidx,drop=FALSE] %*%
            Wmat[ctidx,ct_ref,drop=FALSE]
      }
    }
    
    is_all_zero <- which(unlist(lapply(neighbor_expr_byct, function(x) any(x@x!=0)))==FALSE)
    for(nbz in is_all_zero){
      neighbor_expr_byct[[nbz]] <- .replace_empty_sm(neighbor_expr_byct[[nbz]])
    } 
    rslist <- lapply(neighbor_expr_byct, Matrix::rowSums)
    nz_rows <-  which(unlist(lapply(rslist, function(x) any(x==0))))
    nz_rows <- setdiff(nz_rows, is_all_zero)
    for(nzr in nz_rows){
      neighbor_expr_byct[[nzr]] <- .replace_empty_rows_one_zero(neighbor_expr_byct[[nzr]])
    } 
    cslist <- lapply(neighbor_expr_byct, Matrix::colSums)
    nz_cols <-  which(unlist(lapply(cslist, function(x) any(x==0))))
    for(nzc in nz_cols){
      neighbor_expr_byct[[nzc]] <- .replace_empty_cols_one_zero(neighbor_expr_byct[[nzc]])
    } 
    
    
     
    setnames(metainfo
             ,new=c(cell_type_metadata_colname,cellid_colname)
             ,old=c("cell_type_metadata_colname__"
                    ,"cell_ID"))
   
    if(missing(adjacency_counts_by_ct)){
      ### (possibly weighted) number of adjacent cells by cell type
      adjacency_counts_by_ct <- Matrix::t(Wmat[rownames(ct_matrix),rownames(ct_matrix)]) %*%
        ct_matrix
      cellid <- rownames(adjacency_counts_by_ct) 
      adjacency_counts_by_ct <- as.data.table(as.matrix(adjacency_counts_by_ct))[,cell_ID:=cellid]
    } 
    
    
    return_obj <-  
    list(neighbor_expr_byct=neighbor_expr_byct
                ,adjacency_mat=Wmat
                ,ct_matrix=ct_matrix
                ,adjacency_counts_by_ct = adjacency_counts_by_ct
                ,ref_celltype=ref_celltype
                )
    
    class(return_obj) <- append(class(return_obj), "neighborexpr")
    return(return_obj)
  }


#' Measure gene expression in neighboring cells
#'
#' Measure gene expression in neighboring cells using data.table of cell adjacencies.
#' 
#' @param assay_matrix (normalized) counts matrix.  Normalized is recommended when output of this function is passed to downstream DE to adjust for neighboring cell expression. 
#' @param metadata data.table or data.frame of meta.data associated with assay.  
#' Should contain the cell_type_metadata_colname, and cellid_colname.
#' @param cell_adjacencies data.table of cell connections with columns 'from' and 'to', such as that in Giotto spatial netowork DT object.
#' @param ref_celltype calculate neighbor expression for cells of this particular `ref_celltype`.  
#' @param cell_type_metadata_colname  column name in meta.data of cell types 
#' @param cellid_colname  column name  of cell ids. Or if missing, 
#' @param aggregation "sum" (default), calculates total (weighted) sum amount of neighbor expression.  Otherwise if "average" ,a (weighted) average. 
#' @param weight_colname If non-missing, uses the "weight_colname" column of cell_adjacencies to weight the calculated expression of neighbors. 
#' @param Wmat pre-computed, possibly weighted neighborhood adjacency matrix, calculated based on cell_adjacencies data.table.  
#' If missing, will be computed in this function call.
#' @param ct_matrix pre-computed sparse matrix indicating cell type for each cell.
#' If missing, will be computed in this function call. 
#' @param adjacency_counts_by_ct pre-computed, possibly distance weighted count of neighboring cells by each unique celltype.
#' If missing, will be computed in this function call. 
#' @param nb_celltypes_to_individually_calc vector of cell types for which to individually calculate neighbor expression.
#' default is NULL, which will make a separate list in neighbor_expr_by_ct output with matrix of expression by cell type
#' 
#' @examples
#' library(Seurat)
#' library(data.table); setDTthreads(1)
#' datadir<-system.file("extdata", package="smiDE")
#' sem <- readRDS(paste0(datadir, "/small_nsclc.rds"))
#' metainfo <- data.table::data.table(sem@meta.data)
#'
#' 
#' tiss <- split(metainfo, by="tissue")
#' 
#' ## identify adjacent cells within 0.05 distance of each other within common tissue 
#' cell_adjacency_dt <- 
#'   rbindlist(
#'     lapply(1:length(tiss), function(x){
#'       nb <- fast_make_all_neighbors(tiss[[x]][,.(cell_ID, sdimx, sdimy,cell_type)]
#'                                     ,cellid_col = "cell_ID"
#'                                     ,sdimx_col ="sdimx"
#'                                     ,sdimy_col ="sdimy"
#'                                     ,radius = 0.05) 
#'       return(nb) 
#'     })
#'   )
#' 
#' cell_adjacency_dt[]
#' 
#' ## unweighted, using raw counts
#' nblist_unweighted_counts <- 
#'   measure_neighbor_expression_by_celltype(assay_matrix = sem[["RNA"]]@counts 
#'                                           ,metadata=metainfo
#'                                           ,cell_adjacencies = cell_adjacency_dt
#'                                           ,cell_type_metadata_colname = "cell_type"
#'                                           ,ref_celltype = "mast"
#'                                           ,weight_colname = NULL
#'   )
#' 
#' spp1_neighbor_expr_cts <- 
#'   extract_neighborexpr_by_gene(nblist_unweighted_counts$neighbor_expr_byct, gene= "SPP1")
#' 
#' ## Note that the sum of expression for individual cell types is equal to total 'all cell type' / 'allct' expression.
#' ## Note that the 'otherct' column indicates expression from cell types other than the reference cell type (in this case "mast").
#' spp1_neighbor_expr_cts[1:3]
#' 
#' totalcounts <- sem@meta.data[["totalcounts"]]
#' sem <- Seurat::SetAssayData(sem
#'                            ,"data"
#'                            ,sem[["RNA"]]@counts %*% Matrix::Diagonal(x=mean(totalcounts)/totalcounts, names=colnames(sem))
#'                            )
#' 
#' ## weighted by distance.  Using normalized counts
#' nblist <- 
#'   measure_neighbor_expression_by_celltype(assay_matrix = sem[["RNA"]]@data 
#'                                           ,metadata=metainfo
#'                                           ,cell_adjacencies = cell_adjacency_dt
#'                                           ,cell_type_metadata_colname = "cell_type"
#'                                           ,ref_celltype = "mast"
#'                                           ,weight_colname = "weight"
#'   )
#' 
#' names(nblist)
#' 
#' spp1_neighbor_expr <- 
#'   extract_neighborexpr_by_gene(nblist$neighbor_expr_byct, gene= "SPP1")
#' 
#' spp1_neighbor_expr[1:3]
#' 
#' 
#' 
#' @export
measure_neighbor_expression_by_celltype <- 
  function(assay_matrix
           ,metadata
           ,cell_adjacencies
           ,ref_celltype
           ,cell_type_metadata_colname
           ,aggregation="sum"
           ,cellid_colname="cell_ID"
           ,weight_colname=NULL
           ,Wmat=NULL
           ,ct_matrix = NULL
           ,adjacency_counts_by_ct = NULL
           ,nb_celltypes_to_individually_calc = NULL
  ){
    
    metainfo <- data.table::copy(as.data.table(metadata))
    
    stopifnot(cell_type_metadata_colname %in% names(metainfo)) 
    stopifnot(ref_celltype %in% unique(metainfo[[cell_type_metadata_colname]]))
    
    setnames(metainfo
             , old=c(cell_type_metadata_colname,cellid_colname)
             ,new=c("cell_type_metadata_colname__"
                    ,"cell_ID"))

    
    if(missing(ct_matrix)){
      ct_matrix <- make_ct_matrix(metadata
                                  ,cell_adjacencies[from!=to]
                                  ,cell_type_metadata_colname
                                  ,cellid_colname
      )
    }   
    if(missing(Wmat)){
      Wmat <- make_W(metadata
                     ,cell_adjacencies[from!=to]
                     ,cell_type_metadata_colname
                     ,aggregation
                     ,cellid_colname
                     ,weight_colname
                     )
    } 
      
     
    # extract expression matrix to use, align columns with rows of adjacency matrix for matrix multiplication
    ex <- assay_matrix[,rownames(Wmat),drop=FALSE]
    rm(assay_matrix); gc()
    
    ct_ref <- metainfo[cell_type_metadata_colname__==ref_celltype,cell_ID]
   
    uniq_ct <- metainfo[,unique(cell_type_metadata_colname__)]
    if(is.null(nb_celltypes_to_individually_calc)){
      nb_celltypes_to_individually_calc <- uniq_ct
      nb_celltypes_to_individually_calc <- c("allct", "otherct", nb_celltypes_to_individually_calc) 
    } 
    ind_ct_to_calc <- setdiff(nb_celltypes_to_individually_calc, c("otherct", "allct"))
    stopifnot(all(ind_ct_to_calc %in% uniq_ct))

    ## Neighboring counts by individual cell types 
    neighbor_expr_byct <- vector(mode="list", length = length(nb_celltypes_to_individually_calc))
    names(neighbor_expr_byct) <- gsub("\\-|\\ ","_",nb_celltypes_to_individually_calc)
    
    if("allct" %in% nb_celltypes_to_individually_calc){
      ## All neighboring counts
      neighbor_expr_byct[["allct"]] <- 
        ex %*%
        Wmat[,ct_ref,drop=FALSE]
    }
    
    if("otherct" %in% nb_celltypes_to_individually_calc){
      ## All neighboring counts except those of the same cell type
      neighbor_expr_byct[["otherct"]] <- 
        ex[,setdiff(colnames(ex),ct_ref),drop=FALSE] %*%
        Wmat[setdiff(colnames(ex),ct_ref),ct_ref,drop=FALSE]
      
    } 
    
    if(length(ind_ct_to_calc) > 0){
      for(xx in ind_ct_to_calc){
        ctidx <- metainfo[cell_type_metadata_colname__==xx,cell_ID] 
        neighbor_expr_byct[[gsub("\\-|\\ ","_",xx)]] <-   
          ex[,ctidx,drop=FALSE] %*%
            Wmat[ctidx,ct_ref,drop=FALSE]
      }
    }
    
    is_all_zero <- which(unlist(lapply(neighbor_expr_byct, function(x) any(x@x!=0)))==FALSE)
    for(nbz in is_all_zero){
      neighbor_expr_byct[[nbz]] <- .replace_empty_sm(neighbor_expr_byct[[nbz]])
    } 
    rslist <- lapply(neighbor_expr_byct, Matrix::rowSums)
    nz_rows <-  which(unlist(lapply(rslist, function(x) any(x==0))))
    nz_rows <- setdiff(nz_rows, is_all_zero)
    for(nzr in nz_rows){
      neighbor_expr_byct[[nzr]] <- .replace_empty_rows_one_zero(neighbor_expr_byct[[nzr]])
    } 
    cslist <- lapply(neighbor_expr_byct, Matrix::colSums)
    nz_cols <-  which(unlist(lapply(cslist, function(x) any(x==0))))
    for(nzc in nz_cols){
      neighbor_expr_byct[[nzc]] <- .replace_empty_cols_one_zero(neighbor_expr_byct[[nzc]])
    } 
    
    
     
    setnames(metainfo
             ,new=c(cell_type_metadata_colname,cellid_colname)
             ,old=c("cell_type_metadata_colname__"
                    ,"cell_ID"))
   
    if(missing(adjacency_counts_by_ct)){
      ### (possibly weighted) number of adjacent cells by cell type
      adjacency_counts_by_ct <- Matrix::t(Wmat[rownames(ct_matrix),rownames(ct_matrix)]) %*%
        ct_matrix
      cellid <- rownames(adjacency_counts_by_ct) 
      adjacency_counts_by_ct <- as.data.table(as.matrix(adjacency_counts_by_ct))[,cell_ID:=cellid]
    } 
    
    
    return_obj <-  
    list(neighbor_expr_byct=neighbor_expr_byct
                ,adjacency_mat=Wmat
                ,ct_matrix=ct_matrix
                ,adjacency_counts_by_ct = adjacency_counts_by_ct
                ,ref_celltype=ref_celltype
                )
    
    class(return_obj) <- append(class(return_obj), "neighborexpr")
    return(return_obj)
  }

#' Extract neighboring cell type expression for a particular gene 
#'
#' Extract neighboring cell type expression for a particular gene from an 
#' existing object produced by `smiDE::measure_neighbor_expression_by_celltype`.
#' Columns of this table can be used as covariates in DE, specified in `formula` argument of `run_de` or `smi_de`. 
#' 
#' @param neighbor_expr_list list of matrices by cell type group from `measure_neighbor_expression_by_celltype.`
#' @param gene gene to make data.table for. 
#' 
#' @examples
#' library(Seurat)
#' library(data.table); setDTthreads(1)
#' datadir<-system.file("extdata", package="smiDE")
#' sem <- readRDS(paste0(datadir, "/small_nsclc.rds"))
#' metainfo <- data.table(sem@meta.data)
#' 
#' cell_adjacency_dt <- 
#'   fast_make_all_neighbors(dt = metainfo
#'                           ,cellid_col = "cell_ID"
#'                           ,sdimx_col = "sdimx"
#'                           ,sdimy_col = "sdimy"
#'                           ,radius = 0.05 
#'   )
#' 
#' cell_adjacency_dt[,.(from,to,distance,weight,cell_type_from, cell_type_to)]
#' 
#' totalcounts <- sem@meta.data[["totalcounts"]]
#' sem <- Seurat::SetAssayData(sem
#'                            ,"data"
#'                            ,sem[["RNA"]]@counts %*% Matrix::Diagonal(x=mean(totalcounts)/totalcounts, names=colnames(sem))
#'                            )
#'                            
#' ## weighted by distance.  Using normalized counts
#' nblist <- 
#'   measure_neighbor_expression_by_celltype(assay_matrix =  sem[["RNA"]]@data
#'                                           ,metadata=metainfo
#'                                           ,cell_adjacencies = cell_adjacency_dt
#'                                           ,cell_type_metadata_colname = "cell_type"
#'                                           ,ref_celltype = "mast"
#'                                           ,weight_colname = "weight"
#'   )
#' 
#' spp1_neighbor_expr <- 
#'   extract_neighborexpr_by_gene(nblist$neighbor_expr_byct, gene= "SPP1")
#' 
#' 
#' ## unweighted, using raw counts
#' nblist_unweighted_counts <- 
#'   measure_neighbor_expression_by_celltype(assay_matrix = sem[["RNA"]]@counts
#'                                           ,metadata=metainfo
#'                                           ,cell_adjacencies = cell_adjacency_dt
#'                                           ,cell_type_metadata_colname = "cell_type"
#'                                           ,ref_celltype = "mast"
#'                                           ,weight_colname = NULL
#'   )
#' 
#' spp1_neighbor_expr_cts <- 
#'   extract_neighborexpr_by_gene(nblist_unweighted_counts$neighbor_expr_byct, gene= "SPP1")
#' 
#' ## Note that the sum of expression for individual cell types is equal to total 'all cell type' / 'allct' expression.
#' ## Note that the 'otherct' column indicates expression from cell types other than the reference cell type (in this case "mast").
#' spp1_neighbor_expr_cts[1:3]
#' 
#' @export
#'

extract_neighborexpr_by_gene <- function(neighbor_expr_list,gene, expr_list_names=NULL, cell_IDs=NULL){
  if(is.null(expr_list_names)){
    expr_list_names <- names(neighbor_expr_list)
  }
  if(is.null(cell_IDs)){
    cell_IDs <- colnames(neighbor_expr_list[[1]]) 
  }
  exprl <- 
  lapply(expr_list_names, function(x,cellids=cell_IDs){
     neighbor_expr_list[[x]][gene,cellids]
  })
  names(exprl) <- expr_list_names
  allt <- data.table::data.table(do.call(cbind, exprl))
  return(allt)
}


