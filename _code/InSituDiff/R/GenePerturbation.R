
#' Summarize genes' perturbation levels
#' 
#' For each gene, summarize perturbation with the following statistics:
#' \itemize{
#'  \item How broadly perturbed in disease: its median perturbation in non-controls / its median perturbation in controls
#'  \item How strongly perturbed in disease: its 0.99 quantile perturbation in non-controls / its 0.999 quantile perturbation in controls
#'  \item How variably perturbed in disease: the SD of its perturbation scores in disease
#' }
#' @param x Expression matrix for complete dataset, same as was input to initializeISD. Cells in rows, genes in columns.
#' @param obj Output of initializeISD
#' @param residtype Either "log2ratio" or "diff"
#' @param plotresults Logical, for whether to make a summary plot of the results
#' @param subsetsize Number of cells to use for these calculations
#' @param eps For log2ratio calculations, mean neighborhood expression will be thresholded below at this value
#' @export
summarizeGenePerturbation <- function(x, obj, residtype = "log2ratio", plotresults = FALSE, eps = 1, subsetsize = 1e5) {
  ## format:
  if (is.numeric(obj$tissue)) {
    obj$tissue <- as.character(obj$tissue)
  }
  
  ## checks:
  if (any(is.na(x))) {
    stop("NAs are present in x; complete data is required")
  }
  
  ## get perturbation scores for a subset:
  sub <- InSituDiff:::pseudoRandomSample(vec = seq_len(nrow(x)), n = subsetsize) 
  
  ## calculate perturbation scores:
  mat <- getPerturbations(x = x, obj = obj, cells = sub, genes = NULL, residtype = residtype, eps = eps)
  
  ## summarize genes' perturbation levels:
  sqerr <- mat^2
  
  meanerr <- Matrix::colMeans(sqerr[!obj$iscontrol[sub], ])
  sderr <- sqrt(Matrix::colMeans((sqerr^2)[!obj$iscontrol[sub], ]) - meanerr^2)
  upperconferr <- meanerr + 2 * sderr
  
  meancontrolerr <- Matrix::colMeans(sqerr[obj$iscontrol[sub], ])
  sdcontrolerr <- sqrt(Matrix::colMeans((sqerr^2)[obj$iscontrol[sub], ]) - meancontrolerr^2)
  upperconfcontrolerr <- meancontrolerr + 2 * sdcontrolerr
  
  ratioeps <- 0.1
  out <- cbind((meanerr + ratioeps) / (meancontrolerr + ratioeps), 
               (upperconferr + ratioeps) / (upperconfcontrolerr + ratioeps))
  colnames(out) <- c("broadness", "intensity")
  rownames(out) <- colnames(mat)
  
  if (plotresults) {
    plot(out, col = 0, log = "xy",
         xlab = "Mean disease perturbation / mean control perturbation",
         ylab = "Mean + 2SD disease perturbation / Mean + 2SD control perturbation")
    text(out[, 1], out[, 2], rownames(out), cex = 0.8)
    legend("bottomright", legend = "broadly perturbed", cex = 0.8, bty = "n", text.col = "darkblue")
    legend("topleft", legend = "intensely perturbed", cex = 0.8, bty = "n", text.col = "darkblue")
  }
  
  return(out)
}

#' Identify the most perturbed genes given output of summarizeGenePerturbation()
#' @param genescores Output of summarizeGenePerturbation()
#' @param n Maximum number of genes to return
#' @return A vector of names of highly perturbed genes
#' @export
identifyMostPerturbedGenes <- function(genescores, n = 1000) {
  rownames(genescores)[order(genescores[, "intensity"], decreasing = TRUE)[seq_len(min(n, nrow(genescores) / 5))]]
}


