
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
#' @param eps For log2ratio calculations, this value is added to mean neighborhood expression levels.
#'  Low eps will have higher sensitivity and poorer specificity for detecting perturbations in low expressers. 
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
  modules <- get_modules_from_cor(
    cormat = cormat, 
    min_module_size = 3, max_module_size = 25,          
    resolution = resolution, corthresh = corthresh, min_module_cor = min_module_cor)
  
  return(modules)
}


#' Get coregulation network using soft thresholding
#'
#' Derive a network graph by soft-thresholding the correlation matrix, and computing Topological Overlap Measure (TOM) values
#' @param cormat The correlation matrix
#' @param min_module_size An integer. Won't consider modules smaller than this.
#' @param max_module_size An integer. Not implemented yet.
#' @param resolution Argument to igraph::cluster_leiden. Lower values produce bigger clusters. 
#' @param corthresh Only correlations about this value will go into the adjacency graph fed into leiden clustering
#' @param min_module_cor Only keep modules with average cor above this value.
#' @return An list with two elements. \code{modules}, a list mod module memberships;
#'  and \code{dend}, a dendrogram from hierarchical clustering of the genes
#' @importFrom igraph graph_from_adjacency_matrix
#' @importFrom igraph cluster_leiden
#' @importFrom igraph E
get_modules_from_cor <- function(cormat, min_module_size = 3, max_module_size = 20,
                                 resolution = 0.02, corthresh = 0.1, min_module_cor = 0.1) {
  
  # get adjacency matrix:
  spcor <- cormat
  spcor <- replace(spcor, spcor < corthresh, 0)
  diag(spcor) = 0
  # convert to a graph, with weights = cor^2
  gr <- igraph::graph_from_adjacency_matrix(adjmatrix = spcor^2 * (spcor>0), 
                                            mode = 'undirected', 
                                            weighted = TRUE)
  # leiden cluster:
  leid <- igraph::cluster_leiden(gr, weights = igraph::E(gr)$weight, resolution_parameter = resolution)
  clust <- leid$membership
  names(clust) <- leid$names
  
  # split excessively large clusters - just once:
  for (cid in unique(clust)) {
    genes <- names(clust)[clust == cid]
    if (length(genes) > max_module_size) {
      # subcluster:
      subgraph <- igraph::graph_from_adjacency_matrix(adjmatrix = spcor[genes, genes]^2 * (spcor[genes, genes] > 0), 
                                                      mode = 'undirected', 
                                                      weighted = TRUE)
      subleid <- igraph::cluster_leiden(subgraph, weights = igraph::E(subgraph)$weight, resolution_parameter = resolution * 2)
      subclust <- subleid$membership
      names(subclust) <- subleid$names
      # replace original cluster id:
      clust[genes] = paste0(clust[genes], "subclust", subclust[genes])
    }
  }
  
  # throw out tiny clusters:
  clustersizes <- table(clust)
  clusternames <- names(clustersizes)[clustersizes >= min_module_size]
  
  # throw out clusters with excessively low correlation, and save the rest in a list:
  modules = list()
  for (cid in clusternames) {
    genes = colnames(cormat)[clust == cid]
    meancor = (sum(cormat[genes, genes]) - length(genes)) / (length(genes)^2 - length(genes)) # direct calculation to not break sparse matrix
    
    if (meancor > min_module_cor) {
      newname <- name_module(cormat[genes, genes])
      modules[[newname]] <- genes
    }
  }
  
  modules = modules[order(sapply(modules, length), decreasing = TRUE)]
  return(modules)
}



#' name module based on gene PCs
#'
#' name a module based on its top 2 genes
#' @param mat The conditional correlation matrix for the selected genes
#' @return a name
#' @importFrom Matrix colMeans
name_module <- function(mat) {
  meancors <- Matrix::colMeans(mat)
  n <- ncol(mat)
  top3 <- colnames(mat)[order(meancors, decreasing = TRUE)[1:min(n, 3)]]
  if (n <= 3) {
    name <- paste0(c(top3, n), collapse = "_")
  } else {
    name <- paste0(c(top3[1:2], n), collapse = "_")
  }
  # replace problematic special characters
  name <- make.names(name)
  return(name)
}


