#' Collapse a hierarchy of cell expression profiles to the root level
#' as described in the nested list for cell type hierarchy
#'
#' @param cth Cell type hierarchy, as a nested list with end point cell types
#' in full_profiles
#' @param full_profiles Cell profile matrix, genes x cell types,
#' include all end point cell types in cth
#' @param verbose Option to print progress to console
#' @param print_tree Option to print a structure of the nested lsit
#'
#' @return Matrix of collapsed cell expression profiles
#' @export
#' @import Matrix
#'
#' @examples
#'
#' set.seed(123)
#' test_profiles <- structure(
#'   rnorm(42, mean = 3, sd = 1),
#'   .Dim = c(3L, 14L),
#'   .Dimnames = list(
#'     c("GeneA", "GeneB", "GeneC"),
#'     c(
#'       "B-cell", "endothelial", "fibroblast",
#'       "macrophage", "mast", "mDC", "monocyte", "neutrophil", "NK",
#'       "pDC", "plasmablast", "T4", "T8", "Treg"
#'     )
#'   )
#' )
#'
#' cth_list <- list(
#'   structural = c("endothelial", "fibroblast"),
#'   lymphoid = list(
#'     `B-lymphoid` = c(
#'       "B-cell", "pDC",
#'       "plasmablast"
#'     ),
#'     `T-lymphoid` = c(
#'       "NK", "T4", "T8",
#'       "Treg"
#'     )
#'   ),
#'   myeloid = c(
#'     "macrophage", "mast", "mDC", "monocyte",
#'     "neutrophil"
#'   )
#' )
#'
#' test <- collapseProfiles(cth = cth_list, full_profiles = test_profiles)
#' test
collapseProfiles <- function(cth, full_profiles, verbose = FALSE,
                             print_tree = FALSE) {
  # Verify input match
  if (!all(unlist(cth) %in% colnames(full_profiles))) {
    stop("The 'cth' list has entries not present
         in the colnames of 'full_profiles'")
  }

  # check if there is a list or just the bottom level of the list
  if (class(cth) != "list") {
    return(full_profiles[, cth, drop = FALSE])
  } else {
    if (print_tree) {
      printTree(cth)
    }

    out <- c()

    for (i in seq_along(cth)) {
      if (verbose) {
        message(paste0("Collapsing ", names(cth)[i], " profiles..."))
      }
      if (all(lengths(cth[[i]]) == 1)) {
        # Handle solo endpoints
        if (length(cth[[i]]) == 1 && (names(cth)[i] == "" ||
                                        is.na(names(cth)[i]) ||
                                        is.null(names(cth)[i]))) {
          names(cth)[i] <- cth[[i]]
        }
        pro <- Matrix::rowMeans(full_profiles[, unlist(cth[[i]]), drop = FALSE])
        if (class(out)[1] %in% c("matrix", "data.frame")) {
          out <- cbind(out, matrix(pro,
            ncol = 1,
            dimnames = list(
              rownames(full_profiles),
              names(cth)[i]
            )
          ))
        } else {
          out <- matrix(pro,
            ncol = 1,
            dimnames = list(
              rownames(full_profiles),
              names(cth)[i]
            )
          )
        }
      } else {
        pro <- Matrix::rowMeans(collapseProfiles(
          cth = cth[[i]],
          full_profiles = full_profiles,
          print_tree = FALSE,
          verbose = verbose
        ))
        if (class(out)[1] %in% c("matrix", "data.frame")) {
          out <- cbind(out, matrix(pro,
            ncol = 1, dimnames =
              list(
                rownames(full_profiles),
                names(cth)[i]
              )
          ))
        } else {
          out <- matrix(pro,
            ncol = 1,
            dimnames = list(
              rownames(full_profiles),
              names(cth)[i]
            )
          )
        }
      }
    }

    return(out)
  }
}
