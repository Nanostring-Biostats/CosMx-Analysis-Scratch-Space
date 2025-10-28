
#' add_ncells
#' 
#' @description show the number of cells next to the cluster in a marker heatmap plot.
#' 
#' @export 
add_ncells <- function(heatmap_plot){
  newp <- copy(heatmap_plot)
  lvs <- levels(newp$data$cluster)
  havlvs <- intersect(lvs, unique(as.character(newp$data$cluster)))
  newlvs <- newp$data[,head(.SD,1),by=cluster]
  newlvs <- newlvs[match(havlvs, cluster),paste0(cluster, " (", scales::comma(ncells), ")")]
  newp$data[,cluster:=paste0(cluster, " (", scales::comma(newp$data$ncells), ")")]
  newp$data[,cluster:=factor(cluster, levels=newlvs)]
  return(newp)
}

