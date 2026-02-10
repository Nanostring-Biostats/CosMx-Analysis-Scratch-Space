#' Helper function to demonstrate how to 
#' make a markerslist compatible with `fit_metagene_scores()` function.
#' 
#' Each element of the returned 'markerslist' corresponds to a cell type class, and 
#' Typically, one would want the cell type classes in a 'markerslist' to be mutually exclusive.
#' See examples for more details.
#' 
#' @param index_marker Named list of character vectors specifying the primary marker gene(s) for each cell type class.
#'        Names should correspond to cell type class names.
#' @param predictors Named list of character vectors specifying predictor genes for each cell type class.
#'        Names should match those in \code{index_marker}.
#' @param use_offclass_markers_as_negative_predictors Logical (or named logical vector).
#'        If TRUE, predictors of a cell type class are added as negative predictors
#'        of other celltype classes
#'        (unless they are specified as positive predictors in both classes).
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
    ,use_offclass_markers_as_negative_predictors = TRUE
){
  
  if(is.null(names(index_marker))) stop("index_marker must be a named list, with names corresponding to a celltype class")
  if(is.null(names(predictors))) stop("predictors must be a named list, with names corresponding to a celltype class")
  class_missing_index_marker <- setdiff(names(predictors), names(index_marker))
  class_missing_predictors <- setdiff(names(index_marker), names(predictors))

  stopifnot(length(use_offclass_markers_as_negative_predictors) == 1 || length(use_offclass_markers_as_negative_predictors) == length(predictors))
  
  if(length(use_offclass_markers_as_negative_predictors)==1){
    use_offclass_markers_as_negative_predictors <- rep(use_offclass_markers_as_negative_predictors, length(index_marker))
    names(use_offclass_markers_as_negative_predictors) <- names(predictors)
  } else {
    if(is.null(names(use_offclass_markers_as_negative_predictors))){
      names(use_offclass_markers_as_negative_predictors) <- names(predictors)
    }
  }
  stopifnot(all(names(use_offclass_markers_as_negative_predictors) %in% names(predictors))) 
   
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

  for(ct in names(index_marker)){
    index_markers_not_in_predictors <- setdiff(index_marker[[ct]], predictors[[ct]])
    if(length(index_markers_not_in_predictors) > 0){
      msg <- paste0(paste0(index_markers_not_in_predictors, collapse=","), " specified as index_marker for celltype ", ct
                    ," but not specified as predictors."
                    ,"\nThese will be added to the predictors list and used when not also specified as the response variable in metagene model.")
      message(msg)
      predictors[[ct]] <- unique(c(index_markers_not_in_predictors, predictors[[ct]]))
    }
  }
  
  markerslist <- lapply(names(index_marker), function(ct){
    list(
      index_marker = index_marker[[ct]]
      ,predictors = unique(predictors[[ct]])
      ,use_offclass_markers_as_negative_predictors = use_offclass_markers_as_negative_predictors[[ct]]
    )
  }) 
  names(markerslist) <- names(index_marker)
  class(markerslist) <- append(class(markerslist), "markerslist") 
  return(markerslist)
}


