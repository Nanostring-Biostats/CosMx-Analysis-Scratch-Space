#' Integrate unsupervised and supervised (immune-celltype) cluster labels
#' 
#' @param metadata data.table or data.frame of metadata containing at least c("cell_ID", unsupervised_colname, supervised_colname) columns
#' @param unsupervised_colname colname for the unsupervised (i.e., leiden) clusters
#' @param supervised_colname colname for the supervised (hieratype) clusters
#' @param supervised_labels_keep HieraType cell type categories (i.e., immune celltypes) that should be kept / override the unsupervised labels if they are present.
#' @param adjacency_mat cells x cells sparse similarity matrix used to determine how to dissolve rare clusters 
#' (When a unsupervised cluster is predominantly immune, but HieraType doesnt' label it as such, then there may be leftover cells - 
#' these get relabeled to the most common cluster call amongst that cell's neighbors).  Recommended to use the same adjacency matrix used to generate HieraType results.
#' @param dissolve_smallcluster_if_overwritten_prop_greaterthan This argument is paired with `and_dissolve_smallcluster_if_finalcluster_prop_lessthan`.
#' For an unsupervised cluster which is mostly overwritten (i.e., at least 90% of cells are re-labeled to HieraType), then re-assign the cells in that
#'  cluster to the most common label amongst similar neighboring cells in the adjacency matrix IF the final cluster proportion is sufficiently small.
#'  (i.e., less  than 5% of all cells.).
#'  The defaults below set up two conditions for re-assigning cells in rare unsupervised clusters which were mostly designated to a HieraType cluster.
#'  Example:
#'  `dissolve_smallcluster_if_overwritten_prop_greaterthan = c(0.9, 0.5)` 
#'  `and_dissolve_smallcluster_if_finalcluster_prop_lessthan = c(0.05, 0.01)` 
#'  This means that cells in for particular unsupervised clusters would be re-assigned if
#'  >90% were overwritten by hieratype labels and the final unsupervised cluster proportion < 0.05 or if
#'  >50% were overwritten by hieratype labels and the final unsupervised cluster proportion < 0.01.
#'  As shown in the example, these two arguments should be numeric vectors of equal length.
#' @param and_dissolve_smallcluster_if_finalcluster_prop_lessthan This argument is paired with `dissolve_smallcluster_if_overwritten_prop_greaterthan`.  
#' See notes above for guidance.

#' @return a data.table with ('cellid_colname', 'unsupervised_colname', 'supervised_colname', 'celltype', and 'dissolved_cell'), 
#' where celltype is the new integrated label and dissolved_cell indicates whether the cell came from a dissolved cluster and was re-assigned using nearest neighbors in the adjacency matrix.
celltype_label_integration <- function(metadata
                                       ,adjacency_mat
                                       ,cellid_colname = "cell_ID"
                                       ,unsupervised_colname = "clust"
                                       ,supervised_colname = "celltype_granular"
                                       ,supervised_labels_keep =unique(c(names(HieraType::markerslist_immune)
                                                                         ,names(HieraType::markerslist_tcellmajor)
                                                                         ,names(HieraType::markerslist_cd4tminor)
                                                                         ,names(HieraType::markerslist_cd8tminor))
                                       )
                                       ,dissolve_smallcluster_if_overwritten_prop_greaterthan = c(0.9, 0.5)
                                       ,and_dissolve_smallcluster_if_finalcluster_prop_lessthan = c(0.05, 0.01)
){
  met <- data.table::as.data.table(data.table::copy(metadata))
  stopifnot(unsupervised_colname %in% colnames(met)) 
  stopifnot(supervised_colname %in% colnames(met)) 
  stopifnot(cellid_colname %in% colnames(met)) 
  stopifnot(length(unique(met[[cellid_colname]]))==nrow(met))
 
  if(cellid_colname !="cell_ID" & "cell_ID" %in% colnames(met)) met[["cell_ID"]] <- NULL
  data.table::setnames(met, cellid_colname, "cell_ID")
 
  if(is.null(rownames(adjacency_mat)) | is.null(colnames(adjacency_mat))){
    stop("adjacency_mat needs rownames and colnames corresponding to metadata cell ids.")
  } 
  stopifnot(all(colnames(adjacency_mat) %in% met[["cell_ID"]])) 
  stopifnot(all(rownames(adjacency_mat) %in% met[["cell_ID"]])) 
  
  ## ensure these are aligned 
  if(!all.equal(colnames(adjacency_mat),met[["cell_ID"]])){
    adjacency_mat <- adjacency_mat[,met[["cell_ID"]]] 
  }
  
  if("celltype" %in% colnames(met)){
    data.table::setnames(met, "celltype","celltype_tmp")
    if(supervised_colname=="celltype"){
      supervised_colname <- "celltype_tmp"
    } 
  } 

 
  ## overwrite unsupervised clusters with immune cell labels
  met[["celltype"]] <- data.table::copy(met[[supervised_colname]])
  data.table::set(x = met
                  ,i = met[!celltype %in% supervised_labels_keep,which=TRUE]
                  ,j = "celltype"
                  ,value = met[!celltype %in% supervised_labels_keep][[unsupervised_colname]]
  )
  
  clust <- met[[unsupervised_colname]]
  celltype <- met[["celltype"]]
  names(clust) <- names(celltype) <- met[["cell_ID"]]
  
  propoverwritten <- by(clust != celltype, clust, mean)
  proppercluster <- table(celltype) / length(celltype)
  clusterstodrop <- c()
  for(ii in 1:length(dissolve_smallcluster_if_overwritten_prop_greaterthan)){
    dropii <- 
    intersect(
      names(which(propoverwritten > dissolve_smallcluster_if_overwritten_prop_greaterthan[ii])), ## > 90% of cells were overwritten with a hieratype immune label
      names(which(proppercluster < and_dissolve_smallcluster_if_finalcluster_prop_lessthan[ii])))
    clusterstodrop <- union(clusterstodrop, dropii)
  }
  
  suppressWarnings(
    met[,celltype_idx:=.GRP,by=.(celltype)]
  )
  ct_indicator <- 
    Matrix::sparseMatrix(i=1:nrow(met), j = met$celltype_idx, x = 1)
  dimnames(ct_indicator) <- list(c(met$cell_ID), (met[order(celltype_idx),unique(celltype)]))
  adjacency_mat@x <- rep(1, length(adjacency_mat@x))
  
  nbtab <- adjacency_mat[met[celltype %in% clusterstodrop][["cell_ID"]],,drop=FALSE] %*% ct_indicator
  dimnames(nbtab) <- list(c(met[celltype %in% clusterstodrop][["cell_ID"]]), colnames(ct_indicator))
  nbtab <- as.matrix(nbtab)
  nbtab <- nbtab[, setdiff(colnames(nbtab), clusterstodrop),drop=FALSE]
  relbls <- colnames(nbtab)[apply(nbtab, 1, which.max)]
  celltype[rownames(nbtab)] <- relbls
  met[,celltype:=NULL]
  met[,celltype:=celltype[cell_ID]]
  met[,dissolved_cell:=FALSE]
  met[cell_ID %in% rownames(nbtab),dissolved_cell:=TRUE]
  data.table::setnames(met, "cell_ID", cellid_colname)
  return(met[,c(cellid_colname, unsupervised_colname, supervised_colname, "celltype", "dissolved_cell"),with=FALSE])
}


