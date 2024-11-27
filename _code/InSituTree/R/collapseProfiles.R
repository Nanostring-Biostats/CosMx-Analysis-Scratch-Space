#' Collapse a hierarchy of cell expression profiles according to cell relationship list
#'
#' @param cth Cell type hierarchy, as a nested list
#' @param full_profiles Cell profile matrix, genes x cell types, for all end point cell types
#' @param verbose Option to print progress to console
#' @param print_tree Option to print a structure of the nested lsit
#'
#' @return Matrix of collapsed cell expression profiles
#' @export
#' @import Matrix
#'
#' @examples
#'
#' test_profiles <- structure(c(20.9603386399337, 0.000458563513110102, 0.0213134197201323,
#' 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 13.9034661915883, 0.000495994103622096,
#' 0.0336462476129761, 15.2527213840745, 0.00050482883622113, 23.9519748436983,
#' 47.5121404233421, 0.00107558038855557, 333.727308745673, 13.7004900517181,
#' 0.00037497098661991, 0.00037497098661991, 14.7091889168433, 0.000468241957476138,
#' 0.0288427303176309, 11.9226523626362, 0.00127986493534602, 0.432247094843415,
#' 13.0851834981941, 0.000408071052724823, 0.0245660028716463, 3.1565517220473,
#' 0.000405246861377072, 0.0117737186541241, 13.5470516404443, 0.0003847773659826,
#' 0.126099834470038), 
#' .Dim = c(3L, 14L), 
#' .Dimnames = list(c("FAM138A", "OR4F5", "RNU6-1100P"),
#'     c("B-cell", "endothelial", "fibroblast",
#'     "macrophage", "mast", "mDC", "monocyte", "neutrophil", "NK",
#'     "pDC", "plasmablast", "T4", "T8", "Treg")))
#'
#' cth_list <- list(structural = c("endothelial", "fibroblast")
#'                   , lymphoid = list(`B-lymphoid` = c("B-cell", "pDC", "plasmablast")
#'                                     , `T-lymphoid` = c("NK", "T4", "T8", "Treg"))
#'                   , myeloid = c("macrophage", "mast", "mDC", "monocyte", "neutrophil")
#'                   )
#'
#' test <- collapseProfiles(cth = cth_list, full_profiles = test_profiles)
#' test
collapseProfiles <- function(cth, full_profiles, verbose = F, print_tree = F){
  # Verify input match
  if(!all(unlist(cth) %in% colnames(full_profiles))){
    stop("The 'cth' list has entries not present in the colnames of 'full_profiles'")
  }

  # check if there is a list or just the bottom level of the list
  if(class(cth) != "list"){
    return(full_profiles[, cth, drop = F])
  } else {
    if(print_tree){
      printTree(cth)
    }

    out <- c()

    for(i in 1:length(cth)){
      if(verbose){
        message(paste0("Collapsing ", names(cth)[i], " profiles..."))
      }
      if(all(lengths(cth[[i]]) == 1)){

        # Handle solo endpoints
        if(length(cth[[i]]) == 1 && (names(cth)[i] == "" || is.na(names(cth)[i]) || is.null(names(cth)[i]))){
          names(cth)[i] <- cth[[i]]
        }
        pro <- Matrix::rowMeans(full_profiles[, unlist(cth[[i]]), drop = F])
        if(class(out)[1] %in% c("matrix", "data.frame")){
          out <- cbind(out, matrix(pro, ncol = 1, dimnames = list(rownames(full_profiles), names(cth)[i])))
          } else{
            out <- matrix(pro, ncol = 1, dimnames = list(rownames(full_profiles), names(cth)[i]))
          }
      } else {
        pro <- Matrix::rowMeans(collapseProfiles(cth = cth[[i]], full_profiles = full_profiles, print_tree = F, verbose = verbose))
        if(class(out)[1] %in% c("matrix", "data.frame")){
          out <- cbind(out, matrix(pro, ncol = 1, dimnames = list(rownames(full_profiles), names(cth)[i])))
          } else{
            out <- matrix(pro, ncol = 1, dimnames = list(rownames(full_profiles), names(cth)[i]))
          }
      }
    }

    return(out)
  }

}
