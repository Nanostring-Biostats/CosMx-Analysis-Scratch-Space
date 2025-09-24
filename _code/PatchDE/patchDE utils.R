#' Divide cells into spatially contiguous patches with high statistical power to study a variable of interest.
#' @param xy Matrix of cells' positions
#' @param X Vector to be studied within each district, i.e. vector to optimize power to study.
#' @param dbscan_eps arg for dbscan clustering of hotspot cells
#' @param maxdbsize break up dbclusts bigger than this
#' @param maxradius cells must be at least this close to a hotspot cell to be included in a cluster
#' @return A vector of district assignments. Can be NA. 
drawPatches <- function(xy, X, dbscan_eps = 0.05, maxdbsize = 50, maxradius = 0.3, totvarthresh = 100) {
  
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
  ishot <- gi >= 2
  db <- as.character(dbscan::dbscan(xy[ishot, ], minPts = 1, eps = dbscan_eps)$cluster)
  
  ## subcluster any hotspots of excessive size:
  bigclusts <- names(which(table(db) > maxdbsize))
  for (name in bigclusts) {
    inds <- db == name
    G <- ceiling(sum(inds) / maxdbsize) + 2
    km <- ClusterR::KMeans_rcpp(xy[ishot, ][inds, ], clusters = G, num_init = 1, max_iters = 100)$clusters
    db[inds] <- paste0(db[inds], "_", km)
  }
  
  ## make a neighbors matrix of hotspots based on their locations:
  # get hotspot mean locations:
  hscentroids <- t(sapply(unique(db), function(name) {
    inds <- db == name
    return(c(median(xy[ishot, ][inds, 1]), median(xy[ishot, ][inds, 2])))
  }))
  rownames(hscentroids) <- unique(db)
  hotspotnn <- InSituCor:::nearestNeighborGraph(x = hscentroids[, 1], y = hscentroids[, 2], N = 3) 
  
  ## get the single cells associated with each hotspot:
  closesthot <- FNN::get.knnx(data = xy[ishot, ], 
                              query = xy, 
                              k = 1)
  closesthot$nn.db <- db[closesthot$nn.index]
  closesthot$closeenough <- closesthot$nn.dist < maxradius
  closesthot <- as.data.frame(closesthot)
  
  # for each hotspot, get sum(X) and sum(X2) over the hotspot cell and its connected non-hotspot cells: 
  hs_sumX <- hs_sumX2 <- hs_N <- c()
  for (hs in unique(db)) {
    connectedcells <- which((closesthot$nn.db == hs) & (closesthot$closeenough))
    hs_sumX[hs] <- sum(X[connectedcells])  
    hs_sumX2[hs] <- sum((X^2)[connectedcells])
    hs_N[hs] <- pmax(length(connectedcells), 1)
  }
  
  # run algorithm to agglomerate hotspot groups until power is sufficient
  binseeds <- cluster_by_threshold(adj = hotspotnn, 
                                   X = hs_sumX, 
                                   X2 = hs_sumX2,
                                   N = hs_N,
                                   thresh = totvarthresh)
  
  # map binseeds back to cells in hotspots:
  hsbin <- binseeds[match(db, rownames(hscentroids))]
  
  # and map all cells to a bin:
  cellbin <- hsbin[match(closesthot$nn.db, db)]
  cellbin[!closesthot$closeenough] <- NA
  return(cellbin)
}


# cluster_by_threshold:
#   Inputs:
#     adj    : an n×n sparse symmetric matrix (class "dgCMatrix") of non‐negative edge
#              weights.  adj[i,j] > 0 means i and j are “connected” with weight = adj[i,j].
#     X      : numeric vector length n, sum of X for each node/group
#     X2     : numeric vector length n, sum of X^2 for each node/group
#     N      : numeric vector length n, sum of counts for each node/group
#     thresh : numeric threshold for the group‐metric:
#               metric(C) = sum(X^2)_C − [sum(X)_C]^2 / sum(N)_C
#
#   Behavior:
#     - Start with each node as its own group.
#     - Repeatedly, find all groups with metric < thresh, in ascending order of metric.
#     - For each such group k (smallest metric first), look at all edges (i,j)
#       in adj where i is in group k and j is in a different group m.
#       Pick the edge with smallest weight, and merge group k into group m.
#     - Stop when no group has metric < thresh or no valid merge is possible.
#
#   Returns:
#     Integer vector `cluster_id` of length n, labeling final groups 1..K.
#
cluster_by_threshold <- function(adj, X, X2, N, thresh) {
  if (!inherits(adj, "dgCMatrix")) {
    stop("`adj` must be a sparse dgCMatrix.")
  }
  n <- length(X)
  if (any(lengths(list(X, X2, N)) != n)) {
    stop("Lengths of X, X2, and N must all equal nrow(adj).")
  }
  
  # Precompute the edge list from adj (summary gives i,j, x)
  s <- Matrix::summary(adj)
  # Ensure symmetry: keep both directions or at least one per undirected edge
  # Here we assume adj is symmetric, so summary has both (i,j) and (j,i).
  # We can keep all—merges will consider both.
  
  # Initialize each node in its own cluster
  clusters <- seq_len(n)
  
  # Function to compute metrics for each cluster
  compute_metrics <- function(clusters) {
    ids <- unique(clusters)
    mtrx <- numeric(length(ids))
    names(mtrx) <- ids
    for (id in ids) {
      idx <- which(clusters == id)
      sx  <- sum(X[idx])
      sx2 <- sum(X2[idx])
      sn  <- sum(N[idx])
      mtrx[as.character(id)] <- sx2 - (sx^2) / sn
    }
    mtrx
  }
  
  repeat {
    # 1) Compute current cluster metrics
    metrics <- compute_metrics(clusters)
    # 2) Find clusters below threshold
    below <- names(metrics)[metrics < thresh]
    if (length(below) == 0) break
    
    # 3) Process in ascending order of metric
    below <- below[order(metrics[below])]
    merged_any <- FALSE
    
    for (k_chr in below) {
      k <- as.integer(k_chr)
      # Nodes in cluster k
      idx_k <- which(clusters == k)
      
      # 4) Among edges from idx_k, find those going to a different cluster
      mask <- s$i %in% idx_k & clusters[s$j] != k
      if (!any(mask)) next  # no available neighbor to merge with
      
      # 5) Pick the edge with smallest weight
      cand <- s[mask, , drop = FALSE]
      min_row <- which.min(cand$x)
      j_node  <- cand$j[min_row]
      m       <- clusters[j_node]  # cluster to merge into
      
      # 6) Merge k into m
      clusters[clusters == k] <- m
      merged_any <- TRUE
      break  # restart from top after one merge
    }
    
    # If no merges happened in this pass, we're done
    if (!merged_any) break
  }
  
  # Remap cluster labels to 1..K
  cluster_id <- as.integer(factor(clusters))
  return(paste0("bin", cluster_id))
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


#' patchDE: run DE over all patches
patchDE <- function(y, df, patch) {
  
}