#' Pre-preprocessing for perturbation calculations
#' 
#' Define neighbor and comparator matrices for use in calculating perturbation matrix.
#' Expects a whole study's data, and does fine if passed tissues with diverse disease conditions.
#' If your controls are relevant to different tissues (e.g. healthy liver samples as controls for liver disease
#'  and healthy kidney controls for kidney disease), then run your data in batches 
#'  according to which control tissues are relevant.
#' @param mat Single cell gene expression matrix. Must be non-negative. Normalized counts preferred.
#' @param xy Cells' xy coordinates - a 2-column matrix aligned to the rows of mat
#' @param neighbors Sparse matrix of neighbor relationships. If not provided, xy will be used to calculate this.
#' @param tissue Vector of cells' tissue IDs
#' @param iscontrol Logical vector defining whether each cell belongs to a control tissue / region
#' @param k Number of neighbors to use to define a neighborhood
#' @param radius Radius used to define neighborhoods (choose this or K)
#' @param residtype Either "log2ratio" or "diff". Default to linear-scale "diff".
#' @param makedense Logical, for whether to convert to dense matrix
#' @return A list: \itemize{
#' \item neighbors Matrix giving each cell's spatial neighbors.
#' \item comparators Matrix giving, for each cell, its "pseudoneighbors", i.e.
#'  the cell IDs of its comparator (control) cells' neighbors,
#'  used to concoct a control "pseudoenvironment" for the cell. 
#' \item tissue The tissue vector input to the function.
#' \item iscontrol The iscontrol vector input to the function.
#' }
#' @export
initializeISD <- function(mat, xy, tissue, iscontrol, 
                          k = 100, radius = NULL, neighbor_subsample_rate = 1) {
  
  if (any(mat < 0)) {
    stop("mat must be non-negative")
  }
  if (is.null(rownames(mat))) {
    stop("mat needs rownames")
  }
  ## format:
  if (is.numeric(tissue)) {
    tissue <- as.character(tissue)
  }
  if (sum(iscontrol) == 0) {
    stop("InSituDiff needs some tissues or tissue regions designated as controls")
  }
  # break apart disease and control regions within tissues:
  tissue <- paste0(tissue, "_", c("dis", "ctl")[1 + iscontrol])
  controlnames <- unique(tissue[iscontrol])
  
  # if there's only one control tissue, split it into 5 parts:
  if (length(controlnames) == 1) {
    message("Only 1 control tissue was input. Arbitrarily splitting it into 5 regions. Perturbation values on the control sample will be uncertain; perturbation values of disease samples are unaffected.")
    controlsubsets <- kmeans(xy[iscontrol, ], centers = 5)$cluster
    tissue[iscontrol] <- paste0(tissue[iscontrol], "_region", controlsubsets)
  }
    
  ## checks:
  if (any(is.na(mat))) {
    stop("NAs are present in mat; complete data is required")
  }
  if (any(is.na(xy))) {
    stop("NAs are present in xy; complete data is required")
  }
  if (any(is.na(tissue))) {
    stop("NAs are present in tissue; complete data is required")
  }
  if (any(is.na(iscontrol))) {
    stop("NAs are present in iscontrol; complete data is required")
  }
  
  ## get spatial neighbors for all cells:
  neighbors <- getNeighbors(xy = xy, 
                            neighbors = NULL, 
                            tissue = tissue, 
                            k = k, 
                            radius = radius, 
                            verbose = FALSE) 
  rownames(neighbors) <- rownames(mat)
  colnames(neighbors) <- rownames(mat)
  ## subset neighbors: (not implemented yet - below should be way more efficient)
  if (neighbor_subsample_rate < 1) {
    warning("neighbor_subsample_rate not implemented yet; no subsampling of neighborhoods will be performed")
  }
  
  ## define weights for neighborhood-driven dim reduction: train PCs from neighborhood data subset,
  ## with intent of applying those weights to single cells *before* computing complete neighborhood data.
  withr::with_seed(seed = 0, {sub <- sample(seq_len(nrow(mat)), min(1e5, nrow(mat)))})
  temp_env <- getNeighborhoodExpression(x = mat, 
                                        neighbors = neighbors[sub, ],
                                        makedense = TRUE)   
  # gene scaling factors:
  gene_scaling <- colMeans(temp_env)^-0.5
  # fit PCs:
  wts <- irlba::prcomp_irlba(sweep(temp_env, 2, gene_scaling, "/"), 
                             n = 20, retx = FALSE, scale. = FALSE)$rotation
  # put weights back on the scale of unscaled data:
  wts <- sweep(wts, 1, gene_scaling, "/")
  
  ## get dim-reduced env matrix for all cells:
  neighbormat <- getNeighborhoodExpression(mat %*% wts, 
                                           neighbors = neighbors,
                                           makedense = TRUE)
  
  ## map every cell neighborhood to the closest neighborhood(s) in controls (not of the same sample):
  controlmatches <- matchToControls(x = neighbormat, 
                                    tissue = tissue, 
                                    iscontrol = iscontrol) 
  controlmatches <- matrix(rownames(mat)[controlmatches], nrow(controlmatches))
  rownames(controlmatches) <- rownames(mat)
  
  return(list(neighbors = neighbors,
              controlmatches = controlmatches,
              tissue = tissue,
              iscontrol = iscontrol))
}


#' Calculate perturbations for selected cells and genes
#' 
#' For given cells and genes, compute the matrix of perturbations from matched controls.
#' @param x Expression matrix for complete dataset, same as was input to initializeISD. Cells in rows, genes in columns.
#' @param obj Output of initializeISD
#' @param cells Which cells to produce output for. If NULL, will return results for all cells. 
#' @param genes Which genes to produce output for. To analyze *modules* instead of genes, pass a 
#' named list in which each element holds a vector of gene names for a given module, e.g. 
#' as output by \code{buildGeneModules}. 
#' @param residtype Either "log2ratio" or "diff"
#' @param eps For log2ratio calculations, mean neighborhood expression will be thresholded below at this value.
#'  This is an important parameter: low eps will have higher sensitivity and poorer
#'  specificity for detecting perturbations in low expressers. 
#'  The default of 1 is a somewhat conservative choice. 
#' @return Matrix of perturbation values for the selected cells vs. genes.
#' @export
getPerturbations <- function(x, obj, cells = NULL, genes = NULL, residtype = "log2ratio", eps = 1) {

  checkxandobj(x, obj)
  cells <- formatCellsArg(cells, x)
  controlcells <- as.vector(obj$controlmatches[cells, ])

  # format genes argument:
  genes <- formatGenesArg(genes, x)

  # if a list of metagenes / modules is passed, calculate metagene scores:
  if (is.list(genes)) {
    # calculate module scores
    scores <- c()
    for (i in 1:length(genes)) {
      scores <- cbind(scores, Matrix::rowMeans(x[, genes[[i]]]))
    }
    colnames(scores) <- names(genes)
    # reformat for ingestion by the rest of the function:
    x <- as.matrix(scores)  
    genes <- colnames(x)
    rm(scores)
  }
  
  # warn if too big
  if (length(genes) * length(cells) > 2e9) {
    message(paste0("Use of ", length(cells), " cells and ", length(genes), " genes risks memory issues."))
  }
  
  # compute neighborhood expression for selected genes, cells:
  env <- getNeighborhoodExpression(x[, genes, drop = FALSE], 
                                   neighbors = obj$neighbors[cells, , drop = FALSE],
                                   makedense = TRUE)
  rownames(env) <- cells
  
  # compute neighborhood expression for selected genes, cells:
  controlenv <- getNeighborhoodExpression(x[, genes, drop = FALSE], 
                                   neighbors = obj$neighbors[controlcells, , drop = FALSE],
                                   makedense = TRUE)
  rownames(controlenv) <- controlcells
  
  deltas <- scoreDiff(env = env, 
                      controlenv = controlenv, 
                      residtype = residtype,
                      eps = eps) 
  return(deltas)
}
  

#' Score change from control neighborhoods
#' 
#' Calculates differences between target cell neighborhoods and their matched controls.
#' @param env Environment matrix of target cells
#' @param controlenv Environment matrix of best-matching controls
#' @param residtype Either "log2ratio" or "diff"
#' @return A matrix of changes from nearest controls, in the same dimension as x
#' @export
scoreDiff <- function(env, controlenv, residtype, eps = 1) {
  # get average of matches:
  if (residtype == "log2ratio") {
    return(log2((env + eps) / (controlenv + eps)))     
  }
  if (residtype == "diff") {
    return(env - controlenv)    
  }
}

#' Find the most similar control neighborhood(s) for each cellular neighborhood
#' @param x Matrix of cell x gene neighborhood expression values
#' @param tissue Vector of cells' tissue IDs
#' @param iscontrol Which cells are controls
#' @export
#' @importFrom FNN get.knnx
matchToControls <- function(x, tissue, iscontrol) {
  # build KD tree on controls:
  matches <- matrix(NA, nrow(x), 1)
  # non-control samples compared to nearest control:
  matches[!iscontrol, ] <- which(iscontrol)[
    FNN::get.knnx(data = x[iscontrol, ], query = x[!iscontrol, ], k = 1)$nn.index]  #<--- not sure this will work as intended
  # each control sample compared to all other controls:
  for (name in unique(tissue[iscontrol])) {
    isothercontrol <- iscontrol & (tissue != name) 
    matches[tissue == name, ] <- which(isothercontrol)[
      FNN::get.knnx(data = x[isothercontrol, ], query = x[tissue == name, ], k = 1)$nn.index] 
  }
  return(matches)
}

