#' Tree-style cell typing
#'
#' @param full_profiles Cell profile matrix, genes x cell types, for all end point cell types
#' @param cth Cell type hierarchy, as a nested list
#' @param x Counts matrix, cells x genes
#' @param neg Vector of mean negative controls for each cell
#' @param name_for_new_annotation Name for this annotation
#' @param cohort Vector of cells' cohort membership
#' @param excluded_genes Genes to be excluded during fitting with InSituType
#' @param quantile_absolute_expression_difference_param Quantile cutoff for absolute expression.  Default = 0.5. Both quantile cutoffs must be passed to retain gene.
#' @param quantile_percent_expression_difference_param Quantile cutoff for percent expression.  Default = 0.5. Both quantile cutoffs must be passed to retain gene.
#' @param return_summary_annotation option to return a summary of celltypes.  Default is TRUE.
#'
#' @return List of each InSituType result and a summary of all cell type annotations.
#' @export
#'
#' @examples
#'
#' test_celltype_relationship_list <- list("structural" = 
#'   c("endothelial", "fibroblast")
#'   , "myeloid" = c("macrophage", "mast", "mDC", "monocyte", "neutrophil")
#'   , "lymphoid" = list("B-lymphoid" = c("B-cell", "pDC","plasmablast")
#'     , "T-lymphoid" = list(
#'         "T4" = c("T CD4 memory", "T CD4 naive", "Treg"),
#'         "T8" = c("T CD8 memory", "T CD8 naive"),
#'         "NK")
#'         )
#' )
#'
#' # Extract expression data from insitutype package
#' library(InSituType)
#' data(mini_nsclc)
#' data(ioprofiles)
#'
#'
#' res <- runInSituTree(full_profiles = ioprofiles
#'                    , cth = test_celltype_relationship_list
#'                    , x = as(mini_nsclc$counts, "CsparseMatrix")
#'                    , neg = Matrix::rowMeans(mini_nsclc$neg)
#'                    , name_for_new_annotation = "test"
#'                    , cohort = NULL
#'                    , excluded_genes = c("MALAT1", "B2M", "CD298", 
#'                       "MZT2A", "HLA-A", "HLA-B", "HLA-C")
#'                    )
#'
#'
runInSituTree <- function(x
                          , neg
                          , full_profiles
                          , cth
                          , name_for_new_annotation
                          , cohort = NULL
                          , excluded_genes = c()
                          , quantile_absolute_expression_difference_param = 0.5
                          , quantile_percent_expression_difference_param = 0.5
                          , return_summary_annotation = TRUE
){

  # Argument checks
  if (missing(full_profiles) || missing(cth) || missing(x) || missing(neg) || missing(name_for_new_annotation)) {
    stop("Error: All required arguments (full_profiles, cth, x, neg, name_for_new_annotation) must be provided.")
  }

  if (!is.matrix(full_profiles)) {
    stop("Error: full_profiles must be a matrix.")
  }

  if (!is.list(cth) & !is.vector(cth)) {
    print(paste0("cth is of class ", class(cth)))
    stop("Error: cth must be a list or vector of cell types.")
  }

  if (!is(x, "dgCMatrix")) {
    stop("Error: x must be a sparse matrix of type dgCMatrix.")
  }

  if (!is.numeric(neg)) {
    stop("Error: neg must be a numeric vector.")
  }

  if (!is.character(name_for_new_annotation) || nchar(name_for_new_annotation) == 0) {
    stop("Error: name_for_new_annotation must be a non-empty character string.")
  }

  if (!is.null(cohort) && !is.vector(cohort)) {
    stop("Error: cohort must be a vector if provided.")
  }

  if (!is.null(excluded_genes) && !is.character(excluded_genes)) {
    stop("Error: excluded_genes must be a character vector.")
  }

  if (!is.numeric(quantile_absolute_expression_difference_param) || quantile_absolute_expression_difference_param < 0 || quantile_absolute_expression_difference_param > 1) {
    stop("Error: quantile_absolute_expression_difference_param must be a numeric value between 0 and 1.")
  }

  if (!is.numeric(quantile_percent_expression_difference_param) || quantile_percent_expression_difference_param < 0 || quantile_percent_expression_difference_param > 1) {
    stop("Error: quantile_percent_expression_difference_param must be a numeric value between 0 and 1.")
  }

  out <- list()

  # Make combined profile
  prof <- collapseProfiles(cth = cth, full_profiles = full_profiles)

  # annotate with InSituType
  res <- supervisedSubcluster(reference_profiles = prof,
                              x = x,
                              neg = neg,
                              cohort = cohort,
                              excluded_genes = excluded_genes,
                              quantile_absolute_expression_difference = quantile_absolute_expression_difference_param,
                              quantile_percent_expression_difference = quantile_percent_expression_difference_param
  )
  out[[name_for_new_annotation]] <- list()
  out[[name_for_new_annotation]][["result"]] <- res
  out[[name_for_new_annotation]][["subclusterings"]]<- list()
  out[[name_for_new_annotation]][["name"]] <- name_for_new_annotation

  # annotation next level down

  for(i in names(cth)){
    if(i == "list" || i == name_for_new_annotation || i == ""){
      next}
    message(i)
    celltypes_dropped <- !i %in% unique(res$clust) # check to make sure further subtyping is required.  This is where cell lineages that are not present are dropped out.
    if(celltypes_dropped){
      message(paste0("Warning: No cells found on which to perform further fitting.  Dropping celltypes: ", i))
    }

    if(length(cth[[i]]) > 1 & !celltypes_dropped){
      message(paste0("Annotating  ", i, " cells..."))
      selected_cells <- match( names(res$clust)[res$clust == i], rownames(x) )
      sub_res <- runInSituTree(full_profiles = full_profiles
                               , cth = cth[[i]]
                               , x = x[selected_cells , ]
                               , neg = neg[selected_cells]
                               , name_for_new_annotation = i
                               , cohort = cohort[selected_cells,]
                               , excluded_genes = excluded_genes
                               , quantile_absolute_expression_difference_param = quantile_absolute_expression_difference_param
                               , quantile_percent_expression_difference_param = quantile_percent_expression_difference_param
                               , return_summary_annotation = FALSE
      )
      out[[name_for_new_annotation]][["subclusterings"]] <- append(out[[name_for_new_annotation]][["subclusterings"]], sub_res)
    }
  }

  # return the nested list of insitutype objects
  if(return_summary_annotation){
    summaryAnnotation <- summarizeInSituTree(out[[1]])
    out$summaryAnnotation <- summaryAnnotation

  }

  return(out)

}
