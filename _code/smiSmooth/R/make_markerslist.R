#' Helper function to demonstrate how to 
#' make a markerslist compatible with `fit_metagene_scores()` function.
#' 
#' Each element of the returned 'markerslist' corresponds to a cell type class, and 
#' Typically, one would want the cell type classes in a 'markerslist' to be mutually exclusive.
#' See examples for more details.
#' 
#' @param index_marker 
#' @param predictors 
#' @param signs an optional list of signs for each cell type class, i.e., c(1,1, -1) corresponding to the expected association direction 
#'              of the predictors with the cell type class.  
#'              If left empty (NULL is default), 
#'              all predictors are assumed to be positively associated 
#'              with the cell type class. 
#' @param use_offclass_markers_as_negative_predictors 
#'         if TRUE, predictors of a cell type class are added as negative predictors 
#'         of other celltype classes 
#'         (unless they are specified as positive predictors in both classes).
#'  
#' @return a 'markerslist' list object which can be passed to `fit_metagene_scores`
#'
#' @examples
#'
#' markerslist <-  
#'   make_markerslist(
#'      index_marker = list(bcell = c("CD19")
#'                           ,plasma = c("IGKC")
#'                           ,myeloid = "LYZ"
#'                           ,tcell = "CD3E"
#'                           )
#'      ,predictors = list(tcell = c("CD3D", "CD3E", "CD3G", "CD2", "CD4", "CD8A"  
#'                                   ,"CD8B", "IL2" , "IL4", "IL17A", "IFNG",   "PRF1"  
#'                                   ,"GZMA", "GZMB","FOXP3"
#'                                   )
#'                         ,bcell = c("CD19", "MS4A1", "CD79A", "CD79B", "PAX5"
#'                                      ,"CD22", "IGKC", "TCL1A", "CD40"
#'                                    )
#'                         ,myeloid = c("CD163", "CD68","CSF1R","CD4", "ITGAM", "ITGAX","LYZ"     
#'                                      ,"FCGR3A", "MRC1", "CD14" ,"HLA-DRA", "HLA-DRB1"
#'                                      )
#'                         ,plasma = c("IGHM", "IGHG1", "IGHA1", "IGKC"
#'                                     ,"IGHA1", "IGHD", "CD38")
#'                         )
#'                         
#'    )
#'   
#' @export
#'   
make_markerslist <- function(
    index_marker
    ,predictors
    ,signs = NULL
    ,use_offclass_markers_as_negative_predictors = TRUE
    ,return_as_datatable = FALSE
){
  
  if(is.null(names(index_marker))) stop("index_marker must be a named list, with names corresponding to a celltype class")
  if(is.null(names(predictors))) stop("predictors must be a named list, with names corresponding to a celltype class")
  class_missing_index_marker <- setdiff(names(predictors), names(index_marker))
  class_missing_predictors <- setdiff(names(index_marker), names(predictors))
  if(length(class_missing_index_marker) > 0){
    msg <- paste0("Some of the celltype class names in 'predictors' are not in 'index_marker'"
                  ,"\n"
                  ,paste0(class_missing_index_marker, collapse=", "))
    stop(msg)
  }
  if(length(class_missing_predictors) > 0){
    msg <- paste0("Some of the celltype class names in 'index_marker' are not in 'predictors'"
                  ,"\n"
                  ,paste0(class_missing_predictors, collapse=", "))
    stop(msg)
  }
  
  for(cl in names(predictors)){
    predictors[[cl]] <- setdiff(predictors[[cl]], index_marker[[cl]])
  } 
  
  if(is.null(signs)){
    ### assume all predictors are positive markers for the class
    signs <- lapply(predictors, function(x){
      rep(1, length(x))
    })
    names(signs) <- names(predictors)
  } else {
    stopifnot  
  }
  
  if(use_offclass_markers_as_negative_predictors){
    allmarkers <- c(unlist(predictors), unlist(index_marker))
    for(cl in names(predictors)){
      signs[[cl]] <-  (c(predictors[[cl]]
                         ,setdiff(allmarkers, c(index_marker[[cl]], predictors[[cl]]))) %in% 
                         predictors[[cl]]
      )*2 - 1
      predictors[[cl]] <- c(predictors[[cl]]
                            ,setdiff(allmarkers, c(index_marker[[cl]], predictors[[cl]])))
    }
  }
  
  markerslist <- vector(mode = 'list', length=length(index_marker))
  names(markerslist) <- names(index_marker)
  for(cl in names(index_marker)){
    markerslist[[cl]] <- list(
      index_marker = index_marker[[cl]]
      ,predictors = predictors[[cl]]
      ,signs = signs[[cl]]
    ) 
  }
  
  if(return_as_datatable){
    markerslist <- 
      rbindlist(
        lapply(names(markerslist), function(xx){
          data.table(cell_type = xx
                     ,index_marker = markerslist[[xx]]$index_marker
                     ,predictors = markerslist[[xx]]$predictors
                     ,signs = markerslist[[xx]]$signs
          )
        })
      )
    
  } 
  return(markerslist)
}
