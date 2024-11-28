#' Prints nested list in an easy to read format
#'
#' @param cth Cell type hierarchy, as a nested list
#' @param prefix Prefix for nested list printing, used in recursive call.
#' Recommended to leave alone.
#' @export
#'
#' @examples
#' cth_list <- list(
#'   structural = c("endothelial", "fibroblast"),
#'   lymphoid = list(
#'     `B-lymphoid` = c("B-cell", "pDC", "plasmablast"),
#'     `T-lymphoid` = c("NK", "T4", "T8", "Treg")
#'   ),
#'   myeloid = c("macrophage", "mast", "mDC", "monocyte", "neutrophil")
#' )
#' printTree(cth_list)
#'
printTree <- function(cth, prefix = "") {
  if (is.list(cth)) {
    for (i in seq_along(cth)) {
      if (is.null(names(cth)[i])) {
        cat(prefix, cth[[i]], "\n")
      } else if (names(cth)[i] == "") {
        cat(prefix, cth[[i]], "\n")
      } else {
        cat(prefix, names(cth)[i], "\n")
        printTree(cth[[i]], paste0(prefix, "  "))
      }
    }
  } else {
    for (i in seq_along(cth)) {
      cat(prefix, cth[i], "\n")
    }
  }
}
