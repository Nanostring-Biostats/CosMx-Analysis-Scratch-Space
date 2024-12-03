#' Supervised subclustering using a gene subset and insitutypeML
#'
#' @param reference_profiles Matrix of expression profiles, genes x cell types
#' @param x Counts matrix, cells x genes
#' @param neg Vector of mean negative controls for each cell
#' @param quantile_absolute_expression_difference Minimum absolute expression
#' difference within reference_profiles in terms of quantile level among all
#' genes. Default = 0.5. Both quantile cutoffs must be passed to retain gene.
#' @param quantile_percent_expression_difference Minimum percentage expression
#' difference within reference_profiles in terms of quantile level among all
#' genes. Default = 0.5. Both quantile cutoffs must be passed to retain gene.
#' @param excluded_genes Genes to be excluded during fitting with InSituType
#' @param cohort Vector of cells' cohort membership
#'
#' @import InSituType
#'
#' @return a list with all elements returned by InSituType::insitutypeML() and 
#' additional element `ctsPerCell` for a vector of counts per cell within the 
#' chosen subset of genes.
#' @export

supervisedSubcluster <- function(reference_profiles,
                                 x,
                                 neg,
                                 quantile_absolute_expression_difference = 0.5,
                                 quantile_percent_expression_difference = 0.5,
                                 excluded_genes = NULL,
                                 cohort = NULL) {
  nsprofiles <- reference_profiles[intersect(colnames(x),
                                             rownames(reference_profiles)), ]
  # normalize the profiles
  nsprofiles <- apply(nsprofiles, 2, function(x) x / sum(x) * 100)

  ## Determine which genes have significant differences
  # Absolute difference
  diffs <- sapply(rownames(nsprofiles), function(gene) {
    diff(range(nsprofiles[gene, ]))
  })

  # Percentage difference
  percs <- sapply(rownames(nsprofiles), function(gene) {
    diff(range(nsprofiles[gene, ])) / max(nsprofiles[gene, ])
  })

  # Isolate genes with high absolute and percentage differences
  use_genes <- names(diffs)[diffs >
                              quantile(diffs,
                                       quantile_absolute_expression_difference,
                                       na.rm = TRUE) &
                              percs >
                                quantile(percs,
                                         quantile_percent_expression_difference,
                                         na.rm = TRUE)]

  # Remove genes that should not be considered.
  use_genes <- use_genes[!use_genes %in% excluded_genes]
  message("These are the selected ", length(use_genes),
          " genes for subclustering:")
  dput(use_genes)

  # Warn if number of genes is low
  if (length(use_genes) < 300) {
    warning(paste0(
      "There are fewer than 300 genes being used to subcluster cells into ",
      paste(colnames(reference_profiles), collapse = ", "),
      ". These annotations may be unstable."
    ))
  }

  # verify that all cells have at least 1 count with gene subset
  nonZeroCount_idx <- Matrix::rowSums(x[, use_genes]) >= 1
  if (!all(nonZeroCount_idx)) warning(paste0(sum(!nonZeroCount_idx),
                                             " cells have 0 counts 
                                             in the subclustering gene panel."))
  
  # Convert counts to dgCMatrix
  if(!is(x, "dgCMatrix")){
    x <- convertToDgCMatrix(x)
  }

  # Run supervised InSituType on cell and gene subet
  sup_res <- insitutypeML(x[nonZeroCount_idx, use_genes],
    neg = neg[nonZeroCount_idx],
    reference_profiles = nsprofiles,
    cohort = cohort[nonZeroCount_idx]
  )

  # Add counts per cell with the reduced panel
  sup_res[["ctsPerCell"]] <- rowSums(x[, use_genes])
  names(sup_res$ctsPerCell) <- row.names(x)

  return(sup_res)
}
