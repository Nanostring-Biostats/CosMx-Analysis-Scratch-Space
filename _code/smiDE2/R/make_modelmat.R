
#debugonce(make_modelmat)
#dd <- 
#make_modelmat(obj_ct
#              ,groupVar="niche_custom"
#              ,groupVar_levels=c("immune", "tumor_interface")
#              ,neighborhood_counts = nblist
#              ,simfrom_mf
#              ,target="KRT17"
#              )

make_modelmat <- function(assay_matrix
                          ,metadata
                          ,groupVar
                          ,groupVar_levels
                          ,neighborhood_counts
                          ,formula
                          ,target
                          ,cellid_colname = "cell_ID"
                          ){
  
  metainfo <- copy(as.data.table(metadata))
  
  stopifnot(!is.null(rownames(metadata)) & (cellid_colname %in% names(metainfo)))
  if(!(cellid_colname) %in% names(metainfo)){
    cellid_colname <- "cell_ID"
    metainfo[[cellid_colname]] <- rownames(metadata)
  }
  setnames(metainfo, old=cellid_colname, new="cell_ID")
  
  mTerms <- all.vars(formula)
  mTerms <- setdiff(mTerms, "1")
  neighborTerms <- setdiff(mTerms, names(metainfo))
  candidate_neighborTerms <- unlist(lapply(names(neighborhood_counts$neighbor_expr_byct)
                                           ,function(x) paste0(x, c("", "_intensity", "_expr", "_cellcount"))))
  neighborTerms <- intersect(neighborTerms, candidate_neighborTerms)
  mTerms <- setdiff(mTerms, neighborTerms)
  
  
  
  if(missing(neighborhood_counts)) neighborhood_counts <- list()
  missingTerms <- setdiff(c(mTerms, neighborTerms)
                                  ,c(names(metainfo),names(neighborhood_counts$neighbor_expr_byct))
  )
  
  pDat <- data.table(metainfo)[,c("cell_ID", mTerms),with=F]# sData(object)[,mTerms]
  stopifnot("totalcounts" %in% colnames(metainfo)) 
    
  pDat[,totalcounts:=metainfo$totalcounts]
 
  if(missing(groupVar_levels)){
    groupVar_levels <- unique(pDat[[groupVar]])
  }
  pDat[[groupVar]] <- factor(pDat[[groupVar]], levels=groupVar_levels) 
  for (i in setdiff(names(pDat), groupVar)){
    if (inherits(pDat[[i]], "character")) {
      pDat[, i] <- as.factor(pDat[[i]])
    }
  }
  updatedFormula <- formula(paste("y", as.character(formula)[2]
                                          , sep = " ~ "))
  
  if(!is.null(neighborhood_counts)){
    neighb_dt <- extract_neighborexpr_by_gene(neighborhood_counts$neighbor_expr_byct
                                              ,target)
    
    pDat <- merge(pDat,neighb_dt, by="cell_ID")
    adjacency_ct <- Matrix::t(neighborhood_counts$adjacency_mat[rownames(neighborhood_counts$ct_matrix)
                                                            ,rownames(neighborhood_counts$ct_matrix)]
    ) %*% neighborhood_counts$ct_matrix
    neighbor_ct_count <- data.table(cell_ID=pDat[,cell_ID]
                                    ,as.matrix(adjacency_ct[pDat[,cell_ID],]))
    names(neighbor_ct_count) <- gsub("\\-|\\ ","_",names(neighbor_ct_count))
    
    pDat <- merge(pDat, neighbor_ct_count,suffixes=c("_expr", "_cellcount"))
    all_celltypes <- setdiff(names(neighbor_ct_count), "cell_ID")
    ref_celltype <- neighborhood_counts$ref_celltype 
    other_ctypes <- setdiff(all_celltypes
                            , c(ref_celltype, "allct", "otherct"))
    pDat[,otherct_cellcount:=rowSums(.SD),.SDcols=paste0(other_ctypes, "_cellcount")]
    pDat[,allct_cellcount:=rowSums(.SD),.SDcols=paste0(all_celltypes, "_cellcount")]
    setnames(pDat, old=c("allct", "otherct"), new=c("allct_expr", "otherct_expr"))
  }
   y <- assay_matrix[target,]
   names(y) = colnames(assay_matrix)
   dat <- merge(data.table(cell_ID=names(y), y=y)
                ,pDat, by=c("cell_ID"))
   
   return(dat)
}
