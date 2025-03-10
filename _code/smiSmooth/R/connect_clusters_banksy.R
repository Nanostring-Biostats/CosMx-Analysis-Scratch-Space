


map_to_seed <- function (seed, val){
  max_seed <- max(as.numeric(seed))
  max_val <- max(as.numeric(val))
  max_both <- max(max_seed, max_val)
  con.mat <- table(seed, val)
  cost.mat <- max(con.mat) - con.mat
  matching <- RcppHungarian::HungarianSolver(cost.mat)$pairs
  if (any(matching[, 2, drop = FALSE] == 0)) {
    matching <- matching[matching[, 2, drop = FALSE] != 0, 
    ]
  }
  from <- matching[, 2, drop = FALSE]
  to <- matching[, 1, drop = FALSE]
  unmapped_from <- setdiff(seq(max_val), from)
  unmapped_to <- seq(length(from) + 1, length.out = length(unmapped_from))
  from <- c(from, unmapped_from)
  to <- c(to, unmapped_to)
  val <- factor(val, labels = to[order(from)])
  val <- factor(val, levels = seq_len(max_both))
  table(seed, val)
  val
}



connect_clusters <- function (metadata, clust_nm, map_to = NULL, verbose = TRUE){
  #clust_nm <- clust_nm[grep("^clust", clust_nm)]
  clust_df <- data.table::data.table(metadata)[,clust_nm,with=FALSE]
  clust_df <- clust_df[,lapply(.SD, as.factor)] 
  clust_df <- as.data.frame(clust_df)
  if (is.null(map_to)) {
    n_clust <- apply(clust_df, 2, function(x) length(unique(x)))
    map_order <- order(n_clust)
    for (i in seq(length(map_order) - 1)) {
      curr_many <- map_order[i + 1]
      curr_few <- map_order[i]
      if (verbose) {
        message(clust_nm[curr_many], " --> ", clust_nm[curr_few])
      }
      #clust_df[, curr_many] <- mapToSeed(clust_df[, curr_few], 
      clust_df[, curr_many] <- map_to_seed(clust_df[, curr_few], 
                                         clust_df[, curr_many])
    }
  } else {
    found <- match(map_to, clust_nm)
    if (is.na(found)) 
      stop(map_to, " not found in cluster names /^clust/.")
    seed <- clust_df[, found]
    for (i in setdiff(seq(length(clust_nm)), found)) {
      clust_df[, i] <- mapToSeed(seed, clust_df[, i])
    }
  }
  return(clust_df)
#  metadata[,clust_nm] <- clust_df
  
  #colData(se)[, clust_nm] <- clust_df
  #se
}


