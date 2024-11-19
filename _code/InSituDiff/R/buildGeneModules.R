
#' Find clusters of genes perturbed in similar regions
#' 
#' @param x Expression matrix for complete dataset, same as was input to initializeISD. Cells in rows, genes in columns.
#' @param obj Output of initializeISD
#' @param genes Which genes to produce output for. 
#'  Default of "highlyperturbed" subsets to the most perturbed genes. If NULL, all genes will be used. 
#' @param resolution Resolution parameter for leiden clustering
#' @param corthresh Correlations with absolute value below this will be rounded to zero to save memory
#' @param min_module_cor Modules must have mean correlation of at least this much to be reported
#' @param subsetsize Number of cells to use for these calculations
#' @param residtype Either "log2ratio" or "diff"
#' @param eps For log2ratio calculations, mean neighborhood expression will be thresholded below at this value.
#'Low eps will have higher sensitivity and poorer specificity for detecting perturbations in low expressers. 
#'  The default of 1 is a somewhat conservative choice. 
#' @return A data frame giving module name, gene name, and gene weight, for all genes included in a module.
#' @export
buildGeneModules <- function(x, obj, genes = "highlyperturbed", 
                             resolution = 0.02, corthresh = 0.1, min_module_cor = 0.1,
                             subsetsize = 1e5, residtype = "log2ratio", eps = 1) {
  
  ## handle genes
  if (!is.null(genes)) {
    if (all(genes == "highlyperturbed")) {
      genescores <- summarizeGenePerturbation(x = x, obj = obj, 
                                              residtype = "log2ratio", plotresults = FALSE, 
                                              eps = eps, subsetsize = subsetsize) 
      genes <- identifyMostPerturbedGenes(genescores, n = 1000) 
    }
  } else {
    genes <- colnames(x)
  }
  
  ## get perturbation scores over quasi-random subset:
  sub <- pseudoRandomSample(vec = seq_len(nrow(x)), n = subsetsize) 
  mat <- getPerturbations(x = x, obj = obj, cells = sub, genes = genes, residtype = residtype, eps = eps)
  
  # get correlation matrix of genes' perturbation scores:
  nonzerosd <- which(apply(mat, 2, sd) > 0)
  cormat <- cor(mat[, nonzerosd])
  if (sum(nonzerosd) == 0) {
    stop("All perturbations were 0. Try a lower eps.")
  }
  
  # define modules:
  modules <- InSituCor:::get_modules_from_cor(
    cormat = cormat, 
    min_module_size = 3, max_module_size = 25,          
    resolution = resolution, corthresh = corthresh, min_module_cor = min_module_cor)
  
  return(modules)
}
