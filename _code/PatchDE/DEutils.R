# Functions in this script:
# - hastyDE: very fast DE modeling using matrix algebra and sparse matrices to model all genes at once
# - getPatches: algorithm to divide a tissue into small patches for separate DE runs
# -- cluster_by_threshold: called by drawPatches, not for direct use
# - getPatchPolys: get polygons circling patches, for use in plots
# - patchDE: algorithm to loop hastyDE over patches

# libraries used: Matrix, InSituCor, spdep, FNN, dbscan, ClusterR, alphahull

#' Very fast DE without best practices
#' Runs simple OLS regression on all genes at once using matrix algebra
#' @param y the normalized counts matrix (cells x genes)
#' @param df a data frame of the variables to be modeled
#' @return A list with:
#'   - effect: matrix of effect sizes (genes x predictors)
#'   - se:     matrix of standard errors
#'   - p:      matrix of p-values
#'   - df_resid: residual degrees of freedom
hastyDE <- function(y, df) {
  if (!is.matrix(y) && !inherits(y, "Matrix")) {
    y <- as.matrix(y)
  }
  X <- model.matrix(~ ., df)  # intercept + predictors
  if (nrow(X) != nrow(y)) stop("nrow(y) must equal nrow(df).")
  
  n <- nrow(X); p <- ncol(X)
  coef_names <- colnames(X)[-1]  # drop intercept
  k <- length(coef_names)
  if (k == 0) stop("No predictors found.")
  
  # Crossproducts
  XtX <- Matrix::crossprod(X)
  XtY <- Matrix::crossprod(X, y)
  
  # Inverse (with fallback)
  XtX_inv <- tryCatch(
    chol2inv(chol(XtX)),
    error = function(e) solve(XtX)
  )
  
  # Coefficients: p x G
  B <- XtX_inv %*% XtY
  rownames(B) <- colnames(X)
  
  # Residual variance
  yty   <- Matrix::colSums(y * y)
  RSS   <- yty - Matrix::colSums(B * XtY)
  dfres <- n - p
  sigma2 <- RSS / dfres  # length G
  
  # Standard errors: for each predictor j, sqrt(sigma2 * V_jj)
  Vdiag <- diag(XtX_inv)
  SE <- vapply((2:p), function(j) sqrt(sigma2 * Vdiag[j]), numeric(length(sigma2)))
  colnames(SE) <- coef_names
  
  # Effects and p-values
  Effect <- as.matrix(Matrix::t(B[-1, , drop = FALSE]))   # G x k
  Tstat  <- as.matrix(Effect / SE)
  Pval   <- 2 * pt(abs(Tstat), df = dfres, lower.tail = FALSE)
  
  list(
    effect     = Effect,
    se         = SE,
    p          = Pval,
    df_resid   = dfres
  )
}



#' Assign cells to "patches" for stratified differential expression analysis.
#' Suggested tuning parameters to use: 
#'  1. Control granularity: npatches.
#'  2. Control learning rate: bizesize, then possibly alpha. 
#'  3. Patch size and shape: maxradius, roundness
#'  4. Patch position: initwithhotspots
#' 
#' @param xy Matrix of cells' xy positions
#' @param X Vector giving the values of the design variable. Aligned to the rows (cells) of xy. 
#' @param npatches How many patches to create. 
#' @param bitesize How far a patch can expand its borders in any iteration. Controls the "learning rate"
#' @param maxradius Cells beyond this distance from a patch centroid can't be assigned to it. Used to prevent excessively expansive patches.
#' @param roundness Tuning parameter from 0-1 guiding patches' tendency to have smooth borders. Lower values get better performance (more consistent var(X) per patch) but produce more uneven borders.
#' @param initwithhotspots TRUE or FALSE for whether to initialize patches around var(X) hotspots (TRUE) or just by clustering cells' positions (FALSE). Anecdotally, TRUE performs better.
#' @param n_iters How many iterations to run.
#' @param alpha Number from 0-1. How aggressively patches with low var(X) will grab cells from bigger neighboring patches. 
#' @param effectivezerodist Lower threshold for cells' distance to patch boundaries. 
#' @param plotprogress Logical for whether to show progress over iters
#' @return A vector of cells' patch assigments, possibly including nA's. 
#' @export
#' @importFrom InSituCor nearestNeighborGraph
#' @importFrom InSituCor neighbor_mean
#' @importFrom spdep knearneigh
#' @importFrom spdep knn2nb
#' @importFrom spdep nb2listw
#' @importFrom spdep localG
#' @importFrom stats kmeans
getPatches <- function(xy, X, npatches,
                       bitesize = 0.1,
                       maxradius = 0.5,
                       roundness = 0.5,
                       initwithhotspots = TRUE,
                       n_iters = 25,
                       alpha = 0.5,
                       effectivezerodist = 0.025,
                       plotprogress = FALSE) {
  
  #### preliminaries --------------------------------------
  # scale X:
  X <- X / sd(X, na.rm = TRUE)
  
  #### get initial patch assignments --------------------------------------
  if (initwithhotspots) {
    ## Get local X, X^2, var:
    neighbors <- InSituCor:::nearestNeighborGraph(x = xy[, 1], y = xy[, 2], N = 50)
    cnX <- InSituCor:::neighbor_mean(neighbors = neighbors, x = X)
    cnX2 <- InSituCor:::neighbor_mean(neighbors = neighbors, x = X^2)
    cnvar <- cnX2 - cnX^2
    
    ## get variance hotspot stats:
    knn      <- spdep::knearneigh(xy, k = 8)
    nb       <- spdep::knn2nb(knn)
    # Convert to a spatial‐weights object (row‐standardized)
    lw       <- spdep::nb2listw(nb, style = "W", zero.policy = TRUE)
    # Compute the Getis–Ord G* statistic (Gi*) for each cell
    gi  <- spdep::localG(cnvar, lw, zero.policy = TRUE)
    
    ## call hotspots, and link them with dbscan:
    ishot <- gi >= 1
    # cluster hotspot points to get patch seeds:
    hotspotcluster <- stats::kmeans(xy[ishot, ], centers = npatches, iter.max = 30)$cluster
    seeds <- rep(NA, nrow(xy))
    seeds[ishot] <- paste0("patch", hotspotcluster)
  } else {
    ## simple initial clustering: Mclust on xy alone:
    seeds <- paste0("patch", stats::kmeans(xy, centers = npatches, iter.max = 30)$cluster)
  }

  #### set up iterations: --------------------------------------
  # set up data frame to track cells:
  celldf <- data.frame(X = X)
  rownames(celldf) <- rownames(xy)
  celldf$patch <- rep(NA, nrow(celldf))
  celldf$patch <- seeds

  # intialize patch df:
  uniqueseeds <- setdiff(seeds, NA)
  patchdf <- data.frame(totvar = rep(NA, length(uniqueseeds)),
                        hunger = rep(NA, length(uniqueseeds)))
  rownames(patchdf) = uniqueseeds
  
  # get patch centroids:
  centroids <- c()
  for (name in unique(seeds)) {
    centroids <- rbind(centroids, colMeans(xy[(seeds == name) & !is.na(seeds), ]))
  }
  rownames(centroids) = unique(seeds)
  
  # get each cell's distance to each patch:
  cell2patchproximity <- getproximity2patch(xy, centroids, celldf$patch, bitesize = bitesize, effectivezerodist = effectivezerodist, maxradius = maxradius, roundness = roundness)

  # which cells are just too far to assign:
  celldf$toofar <- Matrix::rowSums(cell2patchproximity != 0) == 0
  
  if (plotprogress) {
    patchcols <- rep(c("#E41A1C","#377EB8","#4DAF4A","#984EA3","#FF7F00",
                       "#FFFF33","#A65628","#F781BF","#999999","#66C2A5",
                       "#FC8D62","#8DA0CB","#E78AC3","#A6D854","#FFD92F",
                       "#E5C494","#B3B3B3","#1B9E77","#D95F02","#7570B3"), 10)[1:nrow(patchdf)]
    names(patchcols) <- rownames(patchdf)
    plot(xy, asp = 1, pch = 16, cex = 0.1, col = "grey80", main = "init seeds")
    points(xy, pch = 16, cex = 0.4, col = patchcols[celldf$patch])
  }
  
  #### iterations: --------------------------------------------------
  for (iter in seq_len(n_iters)) {
    print(iter)
    
    # get patch stats:
    for (name in uniqueseeds) {
      patchdf[name, "totvar"] <- var(celldf$X[celldf$patch == name], na.rm = T) * sum(celldf$patch == name, na.rm = T)
    }
    patchdf[, "hunger"] <- 1 / ((1 - alpha) * mean(patchdf[, "totvar"]) + alpha * patchdf[, "totvar"])
    patchdf$hunger <- patchdf$hunger / sum(patchdf$hunger)
    
    # add cells to nearest (hunger-weighted) patch:
    celldf$patch <- rownames(patchdf)[
      apply(cell2patchproximity[, rownames(patchdf)] %*% diag(patchdf$hunger), 1, which.max)]
    celldf$patch[celldf$toofar] = NA
    
    # update patch centroids:
    centroids <- c()
    patchnames <- setdiff(unique(celldf$patch), NA)
    for (name in patchnames) {
      centroids <- rbind(centroids, colMeans(xy[(celldf$patch == name) & !is.na(celldf$patch), ]))
    }
    rownames(centroids) = patchnames
    
    # update distances:
    cell2patchproximity <- getproximity2patch(xy, centroids, 
                                              patch = celldf$patch, 
                                              bitesize = bitesize, 
                                              effectivezerodist, 
                                              maxradius = maxradius,
                                              roundness = roundness)
    celldf$toofar <- Matrix::rowSums(cell2patchproximity != 0) == 0
    
    # plot:
    if (plotprogress) {
      par(mfrow = c(2,1))
      plot(xy, asp = 1, pch = 16, cex = 0.1, col = "grey80", main = iter)
      points(xy, pch = 16, cex = 0.4, col = patchcols[celldf$patch])
      
      barplot(patchdf$totvar, col = patchcols[rownames(patchdf)], ylim = c(0,1000), ylab = "Predictor sum of squares", xlab = "Patches")
    }
  }
  out <- celldf$patch
  names(out) <- rownames(celldf)
  return(out)
}



#' Score every cell's proximity to every patch. 
#' "Proximity" is roughly an inverse distance. Cells sufficiently far from a patch get proximity = 0.
#'  Distance is calculated as a weighted geometric mean of a cell's distance to a patch's center
#'   and to its nearest cell within the patch (i.e. to its border).
#' @param xy Matrix of cells' xy positions
#' @param centroids Matrix of patch centroids (patches * xy position)
#' @param patch Vector of cells' patch assignments. Aligned to the rows to xy.
#' @param bitesize How far a patch can expand its borders in any iteration. Controls the "learning rate"
#' @param maxradius Cells beyond this distance from a patch centroid can't be assigned to it. Used to prevent excessively expansive patches.
#' @param roundness Tuning parameter from 0-1 guiding patches' tendency to have smooth borders. Lower values get better performance (more consistent var(X) per patch) but produce more uneven borders.
#' @param effectivezerodist Lower threshold for cells' distance to patch boundaries. 
#' @return A matrix of cell * patch proximity scores, usually zero for cells that are too far from a given patch.
getproximity2patch <- function(xy, centroids, patch, 
                               bitesize, maxradius, roundness, effectivezerodist) {
  patchnames <- setdiff(unique(patch), NA)
  roundness <- max(min(roundness, 1), 0)
  
  # dist to centroid:  
  centroiddistmat <- matrix(NA, nrow(xy), length(patchnames),
                            dimnames = list(rownames(xy), patchnames))
  for (name in patchnames) {
    centroiddistmat[, name] <- sqrt((xy[, 1] - centroids[name, 1])^2 + (xy[, 2] - centroids[name, 2])^2)
  }
  
  # dist to any cells in patch:
  anydistmat <- matrix(NA, nrow(xy), length(patchnames),
                       dimnames = list(rownames(xy), patchnames))
  for (name in patchnames) {
    anydistmat[, name] <- FNN::get.knnx(
      data  = as.matrix(xy[(patch == name) & !is.na(patch), , drop=FALSE]),
      query = as.matrix(xy),
      k = 1
    )$nn.dist[, 1]
  }
  anydistmat[anydistmat < effectivezerodist] <- effectivezerodist
  
  # weighted average (geomean) of distances to border and centroid:
  bigdist <- exp((1 - roundness) * log2(anydistmat) + roundness * log2(centroiddistmat[rownames(anydistmat), colnames(anydistmat)]))

  # convert to proximity: 
  proxmat <- 1 / bigdist
  
  # censor cells that are too far:
  proxmat[anydistmat > bitesize] <- 0
  proxmat[centroiddistmat[rownames(anydistmat), colnames(anydistmat)] > maxradius] <- 0
  
  return(proxmat)
}



#' Find rough polygon boundaries of patches for visualizations
#' @param xy Cells' xy positions
#' @param patch Vector of patch assignments, aligned to the rows of xy
#' @return A named list of polygons, one per patch
#' @export
getPatchPolys <- function(xy, patch) {
  polys <- list()
  cluster_levels <- unique(patch)
  for (i in seq_along(unique(patch))) {
    k <- cluster_levels[i]
    idx <- which(patch == k)
    
    # Only attempt hull if >= 3 points
    if (length(idx) >= 3) {
      pts_k <- xy[idx, , drop = FALSE]       # M_k × 2 matrix of points in cluster k
      hull_indices <- chull(pts_k)           # indices (1..M_k) along the convex hull
      polys[[i]]  <- pts_k[hull_indices, ]       # hull vertices, in order
      names(polys)[i] <- unique(patch)[i]
    }
  }
  return(polys)
}



#' Get polygon borders of patches for visualizations
#' Uses alphahull::ashape()
#' @param xy Cells' xy positions
#' @param patch Vector of patch assignments, aligned to the rows of xy
#' @return A named list of alphahull::ashape objects, one per patch
#' @export
getPatchHulls <- function(xy, patch, alpha = 0.1) {
  hulls <- list()
  cluster_levels <- unique(patch)
  for (i in seq_along(unique(patch))) {
    k <- cluster_levels[i]
    idx <- which(patch == k)
    
    # Only attempt hull if >= 3 points
    if (length(idx) >= 3) {
      pts_k <- xy[idx, , drop = FALSE]       # M_k × 2 matrix of points in cluster k
      hulls[[i]] <- alphahull::ashape(pts_k, alpha = alpha)
      names(hulls)[i] <- unique(patch)[i]
    }
  }
  return(hulls[sapply(hulls, length) > 0])
}

#' patchDE: run DE over all patches
#' @param y Expression matrix, cells * genes
#' @param df Data frame to be used as DE predictors
#' @param patch Vector of patch IDs
#' @export
patchDE <- function(y, df, patch) {
  # get DE results per patch:
  results <- list()
  patches <- setdiff(unique(patch), NA)
  for (patchid in patches) {
    patchinds <- (patch == patchid) & !is.na(patch)
    results[[patchid]] <- hastyDE(y = y[patchinds, , drop = FALSE],
                                  df = df[patchinds, , drop = FALSE])
  }
  # reformat to a per-variable list:
  variables <- colnames(results[[1]][[1]])
  out <- list()
  for (varname in variables) {
    out[[varname]] <- list()
    out[[varname]]$pvals <- sapply(results, function(tmp){tmp$p})
    out[[varname]]$ests <- sapply(results, function(tmp){tmp$effect})
    out[[varname]]$ses <- sapply(results, function(tmp){tmp$se})
    rownames(out[[varname]]$pvals) <- rownames(out[[varname]]$ests) <- rownames(out[[varname]]$ses) <- rownames(results[[1]][[1]])
  }
  return(out)
}


#' #' Divide cells into spatially contiguous patches with high statistical power to study a variable of interest.
#' #' @param xy Matrix of cells' positions
#' #' @param X Vector to be studied within each district, i.e. vector to optimize power to study.
#' #' @param dbscan_eps arg for dbscan clustering of hotspot cells
#' #' @param maxdbsize break up dbclusts bigger than this
#' #' @param maxradius cells must be at least this close to a hotspot cell to be included in a cluster
#' #' @return A vector of district assignments. Can be NA. 
#' #' @export
#' getPatchesOLD <- function(xy, X, dbscan_eps = 0.05, maxdbsize = 50, maxradius = 0.3, totvarthresh = 100) {
#'   
#'   ## scale X:
#'   X <- X / sd(X, na.rm = TRUE)
#'   
#'   ## Get local X, X^2, var:
#'   neighbors <- InSituCor:::nearestNeighborGraph(x = xy[, 1], y = xy[, 2], N = 50)
#'   cnX <- InSituCor:::neighbor_mean(neighbors = neighbors, x = X)
#'   cnX2 <- InSituCor:::neighbor_mean(neighbors = neighbors, x = X^2)
#'   cnvar <- cnX2 - cnX^2
#'   
#'   ## get variance hotspot stats:
#'   knn      <- spdep::knearneigh(xy, k = 8)
#'   nb       <- spdep::knn2nb(knn)
#'   # Convert to a spatial‐weights object (row‐standardized)
#'   lw       <- spdep::nb2listw(nb, style = "W", zero.policy = TRUE)
#'   # Compute the Getis–Ord G* statistic (Gi*) for each cell
#'   gi  <- spdep::localG(cnvar, lw, zero.policy = TRUE)
#'   
#'   ## call hotspots, and link them with dbscan:
#'   ishot <- gi >= 2
#'   db <- as.character(dbscan::dbscan(xy[ishot, ], minPts = 1, eps = dbscan_eps)$cluster)
#'   
#'   ## subcluster any hotspots of excessive size:
#'   bigclusts <- names(which(table(db) > maxdbsize))
#'   for (name in bigclusts) {
#'     inds <- db == name
#'     G <- ceiling(sum(inds) / maxdbsize) + 2
#'     km <- ClusterR::KMeans_rcpp(xy[ishot, ][inds, ], clusters = G, num_init = 1, max_iters = 100)$clusters
#'     db[inds] <- paste0(db[inds], "_", km)
#'   }
#'   
#'   ## make a neighbors matrix of hotspots based on their locations:
#'   # get hotspot mean locations:
#'   hscentroids <- t(sapply(unique(db), function(name) {
#'     inds <- db == name
#'     return(c(median(xy[ishot, ][inds, 1]), median(xy[ishot, ][inds, 2])))
#'   }))
#'   rownames(hscentroids) <- unique(db)
#'   hotspotnn <- InSituCor:::nearestNeighborGraph(x = hscentroids[, 1], y = hscentroids[, 2], N = 3) 
#'   # sever when distances > radius:
#'   hotspotnn[hotspotnn < maxradius] <- 0
#'   
#'   ## get the single cells associated with each hotspot:
#'   closesthot <- FNN::get.knnx(data = xy[ishot, ], 
#'                               query = xy, 
#'                               k = 1)
#'   closesthot$nn.db <- db[closesthot$nn.index]
#'   closesthot$closeenough <- closesthot$nn.dist < maxradius
#'   closesthot <- as.data.frame(closesthot)
#'   
#'   # for each hotspot, get sum(X) and sum(X2) over the hotspot cell and its connected non-hotspot cells: 
#'   hs_sumX <- hs_sumX2 <- hs_N <- c()
#'   for (hs in unique(db)) {
#'     connectedcells <- which((closesthot$nn.db == hs) & (closesthot$closeenough))
#'     hs_sumX[hs] <- sum(X[connectedcells])  
#'     hs_sumX2[hs] <- sum((X^2)[connectedcells])
#'     hs_N[hs] <- pmax(length(connectedcells), 1)
#'   }
#'   
#'   # run algorithm to agglomerate hotspot groups until power is sufficient
#'   binseeds <- cluster_by_threshold(adj = hotspotnn, 
#'                                    X = hs_sumX, 
#'                                    X2 = hs_sumX2,
#'                                    N = hs_N,
#'                                    thresh = totvarthresh)
#'   
#'   # map binseeds back to cells in hotspots:
#'   hsbin <- binseeds[match(db, rownames(hscentroids))]
#'   
#'   # and map all cells to a bin:
#'   cellbin <- hsbin[match(closesthot$nn.db, db)]
#'   cellbin[!closesthot$closeenough] <- NA
#'   return(cellbin)
#' }


# # cluster_by_threshold:
# #   Inputs:
# #     adj    : an n×n sparse symmetric matrix (class "dgCMatrix") of non‐negative edge
# #              weights.  adj[i,j] > 0 means i and j are “connected” with weight = adj[i,j].
# #     X      : numeric vector length n, sum of X for each node/group
# #     X2     : numeric vector length n, sum of X^2 for each node/group
# #     N      : numeric vector length n, sum of counts for each node/group
# #     thresh : numeric threshold for the group‐metric:
# #               metric(C) = sum(X^2)_C − [sum(X)_C]^2 / sum(N)_C
# #
# #   Behavior:
# #     - Start with each node as its own group.
# #     - Repeatedly, find all groups with metric < thresh, in ascending order of metric.
# #     - For each such group k (smallest metric first), look at all edges (i,j)
# #       in adj where i is in group k and j is in a different group m.
# #       Pick the edge with smallest weight, and merge group k into group m.
# #     - Stop when no group has metric < thresh or no valid merge is possible.
# #
# #   Returns:
# #     Integer vector `cluster_id` of length n, labeling final groups 1..K.
# #
# cluster_by_threshold <- function(adj, X, X2, N, thresh) {
#   if (!inherits(adj, "dgCMatrix")) {
#     stop("`adj` must be a sparse dgCMatrix.")
#   }
#   n <- length(X)
#   if (any(lengths(list(X, X2, N)) != n)) {
#     stop("Lengths of X, X2, and N must all equal nrow(adj).")
#   }
#   
#   # Precompute the edge list from adj (summary gives i,j, x)
#   s <- Matrix::summary(adj)
#   # Ensure symmetry: keep both directions or at least one per undirected edge
#   # Here we assume adj is symmetric, so summary has both (i,j) and (j,i).
#   # We can keep all—merges will consider both.
#   
#   # Initialize each node in its own cluster
#   clusters <- seq_len(n)
#   
#   # Function to compute metrics for each cluster
#   compute_metrics <- function(clusters) {
#     ids <- unique(clusters)
#     mtrx <- numeric(length(ids))
#     names(mtrx) <- ids
#     for (id in ids) {
#       idx <- which(clusters == id)
#       sx  <- sum(X[idx])
#       sx2 <- sum(X2[idx])
#       sn  <- sum(N[idx])
#       mtrx[as.character(id)] <- sx2 - (sx^2) / sn
#     }
#     mtrx
#   }
#   
#   repeat {
#     # 1) Compute current cluster metrics
#     metrics <- compute_metrics(clusters)
#     # 2) Find clusters below threshold
#     below <- names(metrics)[metrics < thresh]
#     if (length(below) == 0) break
#     
#     # 3) Process in ascending order of metric
#     below <- below[order(metrics[below])]
#     merged_any <- FALSE
#     
#     for (k_chr in below) {
#       k <- as.integer(k_chr)
#       # Nodes in cluster k
#       idx_k <- which(clusters == k)
#       
#       # 4) Among edges from idx_k, find those going to a different cluster
#       mask <- s$i %in% idx_k & clusters[s$j] != k
#       if (!any(mask)) next  # no available neighbor to merge with
#       
#       # 5) Pick the edge with smallest weight
#       cand <- s[mask, , drop = FALSE]
#       min_row <- which.min(cand$x)
#       j_node  <- cand$j[min_row]
#       m       <- clusters[j_node]  # cluster to merge into
#       
#       # 6) Merge k into m
#       clusters[clusters == k] <- m
#       merged_any <- TRUE
#       break  # restart from top after one merge
#     }
#     
#     # If no merges happened in this pass, we're done
#     if (!merged_any) break
#   }
#   
#   # Remap cluster labels to 1..K
#   cluster_id <- as.integer(factor(clusters))
#   return(paste0("bin", cluster_id))
# }


