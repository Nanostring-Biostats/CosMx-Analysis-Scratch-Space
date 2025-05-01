#' Cluster cells by their environment perturbations
#' 
#' Finds clusters of neighborhoods with similar perturbation patterns. 
#' 
#' @param x Expression matrix for complete dataset, same as was input to initializeISD. Cells in rows, genes in columns.
#' @param obj Output of initializeISD
#' @param genes Which genes to produce output for. 
#'  Default of "highlyperturbed" subsets to the most perturbed genes. If NULL, all genes will be used. 
#' @param cells Which cells to produce output for. If NULL, will return results for all cells. 
#' @param nclust How many clusters to fit
#' @param residtype Either "log2ratio" or "diff"
#' @param eps For log2ratio calculations, mean neighborhood expression will be thresholded below at this value
#' @importFrom mclust Mclust
#' @importFrom mclust mclustBIC
#' @export
clusterPerturbations <- function(x, obj, cells = NULL, genes = "highlyperturbed", nclust = 12, residtype = "log2ratio", eps = 1) {
  
  ## handle genes
  if (!is.null(genes)) {
    if (all(genes == "highlyperturbed")) {
      genescores <- summarizeGenePerturbation(x = x, obj = obj, 
                                              residtype = "log2ratio", plotresults = FALSE, 
                                              eps = eps, subsetsize = 2000) 
      genes <- identifyMostPerturbedGenes(genescores, n = 1000) 
    }
  } else {
    genes <- colnames(x)
  }
  
  # checks:
  checkxandobj(x, obj)
  cells <- formatCellsArg(cells, x)
  genes <- formatGenesArg(genes, x)
  
  # calculate perturbations:
  perturbations <- getPerturbations(x = x, obj = obj, 
                                    cells = cells, genes = genes, 
                                    residtype = residtype,
                                    eps = eps) 
  # remove 0-variance genes:
  nonzerosd <- which(apply(perturbations, 2, sd) > 0)
  perturbations <- perturbations[, nonzerosd]
  if (sum(nonzerosd) == 0) {
    stop("All perturbations were 0. Try a lower eps.")
  }

  # cluster them: (use leiden on a subset, then take centroids from it)
  mc <- mclust::Mclust(data = perturbations, G = nclust, modelNames = "EII")
  means <- mc$parameters$mean
  rownames(means) <- colnames(perturbations)
  colnames(means) <- paste0("domain", seq_len(ncol(means)))
  return(list(clust = paste0("domain", mc$classification), means = means))
}




