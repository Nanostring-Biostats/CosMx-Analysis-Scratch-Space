#' Find node names of a nested list
#'
#' @param cth Cell type hierarchy, as a nested list
#' @return vector of node names
#' @export
#' @examples
#'
#' cth_list <- list(structural = c("endothelial", "fibroblast")
#'                   , lymphoid = list(`B-lymphoid` = c("B-cell", "pDC", "plasmablast")
#'                                     , `T-lymphoid` = c("NK", "T4", "T8", "Treg"))
#'                   , myeloid = c("macrophage", "mast", "mDC", "monocyte", "neutrophil")
#'                   )
#' findTreeNodeNames(cth_list)
#'
#'
findTreeNodeNames <- function(cth){
  if("list" %in% class(cth)){
    ns <- names(cth)
    out <- c()
    for(i in ns) out <- append(out, c(i, unlist(findTreeNodeNames(cth[[i]]))))
  } else {
    return(NULL)
  }
  return(out)
}
