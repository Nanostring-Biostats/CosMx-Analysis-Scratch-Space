
#' K-means cluster in x-y space, used for downstream spatial random effect models.
#'
#' @param metadata metadata which includes cell-level x and y coordinates 
#' with names corresponding to `x_coord_col` and `y_coord_col` arguments.
#' @param x_coord_col name of column with x-coordinate for cell, i.e., "sdimx".
#' @param y_coord_col name of column with x-coordinate for cell, i.e., "sdimy".
#' @param cluster_suffix output x,y cluster centroids with this suffix, 
#' i.e., "_cluster" => "sdimx_cluster", "sdimy_cluster"
#' @param cluster_name name for column to indicate each unique k-means cluster.
#' @param k  number of kmeans clusters
#' @param k_prop_n  number of kmeans clusters as a function of total number of cells.
#' i.e., with k_prop_n = 0.2 and 100 total cells, return 20 clusters.
#' @param seed seed set for ensuring reproducible k-means results.  
#' .Random.seed state is restored within the function, so global environment is not affected. 
#' @param split_neighbors_by_colname In the case of multiple flow-cells, with overlapping range of x and y-coordinates, 
#' specifying this column ensures that different samples/flow-cells with similar x/y coordinates will not belong to the same cluster.
#'
#' 
#' @export 
xy_kmeans_clusters <- function(metadata
                               ,x_coord_col = "sdimx"
                               ,y_coord_col = "sdimy"
                               ,cluster_suffix = "_cluster"
                               ,cluster_name = "k_cluster"
                               ,k=NULL
                               ,k_prop_n = 0.25
                               ,seed = 1
                               ,split_neighbors_by_colname = NULL
){
  
  met <- copy(as.data.table(metadata))
  stopifnot(all(c(x_coord_col, y_coord_col) %in% colnames(met)))
  if(!missing(split_neighbors_by_colname)){
    stopifnot(split_neighbors_by_colname %in% colnames(met))
  }
  clusterx_name <- paste0(x_coord_col, cluster_suffix)
  clustery_name <- paste0(y_coord_col, cluster_suffix)
  if(clusterx_name %in% colnames(met)){
    warning(paste0(clusterx_name, " already a column in metadata. Overwriting existing column with new clusters."))
  }
  if(clustery_name %in% colnames(met)){
    warning(paste0(clustery_name, " already a column in metadata. Overwriting existing column with new clusters."))
  }
  
  if(!is.null(k) & !(is.null(k_prop_n) | is.na(k_prop_n))){
    msg <- paste0("Warning: using provided number of clusters 'k'=",k, " and ignoring k_prop_n=", k_prop_n)
    warning(msg) 
  }
  if(is.null(k) & is.null(k_prop_n)){
    stop("Either 'k' or 'k_prop_n' argument must be provided.")
  }
  
  if(is.null(k)){
    k <- floor(nrow(met)*k_prop_n)
  }   
  for(nm in c(clusterx_name, clustery_name, cluster_name)){
    if(nm %in% colnames(met)) met[[nm]] <- NULL   
  }
  
  old_seed <- .Random.seed    
  on.exit({ .Random.seed <<- old_seed }) 
  set.seed(seed)
  
  nogroup <- FALSE 
  if(is.null(split_neighbors_by_colname)){
    split_neighbors_by_colname <- "splitcol__"
    met[[split_neighbors_by_colname]] <- "grp"
    nogroup <- TRUE
  } 
  met <- 
    rbindlist(
      lapply(split(met, by=c(split_neighbors_by_colname))
             ,function(xx, kk=k){
               kx <- round(nrow(xx)/nrow(met)*kk)
               km <- kmeans(xx[,c(x_coord_col, y_coord_col),with=FALSE], centers=kx)
               kmclus <- data.table(km$centers)
               kmclus[[cluster_name]] <- 1:nrow(kmclus)
               setnames(kmclus, c(x_coord_col, y_coord_col), c(clusterx_name, clustery_name))
               xx[[cluster_name]] <- km$cluster
               xx <- merge(xx, kmclus, by=cluster_name, sort=FALSE)
               return(xx)
             })  
    ) 
  met[,(c(cluster_name)) := .GRP,by=c(split_neighbors_by_colname, clusterx_name, clustery_name, cluster_name)]
  if(nogroup){
    met[[split_neighbors_by_colname]] <- NULL
  }
  return(met) 
}

