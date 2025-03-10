
#' Make a barplot for proportion of cells expressing marker genes at 
#' different confidence thresholds.
#' 
#' @param normed expression data used to make barplot
#' @param metadata data.frame or data.table containing `score_column`, `cellid_column`, `cluster_column`
#' @param clusters clusters to include in barplots
#' @param markers genes to include in barplots
#' @param thresholds confidence thresholds to evaluate
#' @param score_column column in `metadata` containing the confidence scores.
#' @param cluster_column column in `metadata` containing the clusters.
#' @param cluster_order optional order for the clusters in facet.
#' @param marker_order optional order for the markers to be plotted.
#' 
#' @export
marker_threshold_plot <- function(normed=NULL
                                  ,metadata
                                  ,clusters = NULL
                                  ,markers = NULL
                                  ,thresholds = c(0.5, 0.7, 0.8, 0.9,0.99, 0.999)
                                  ,score_column = "best_score"
                                  ,cellid_column = "cell_ID"
                                  ,cluster_column = "celltype"
                                  ,cluster_order = NULL
                                  ,marker_order = NULL
                                  ){

  if(is.null(markers)) markers <- rownames(normed)
  marker_order <- markers 
  
  met <- data.table::copy(data.table::data.table(metadata))
  if(score_column!="best_score" & "best_score" %in% names(met)){
    met[["best_score"]] <- NULL
  } 
  data.table::setnames(met, score_column, "best_score")
  
  if(cellid_column!="cell_ID" & "cell_ID" %in% names(met)){
    met[["cell_ID"]] <- NULL
  } 
  data.table::setnames(met, cellid_column, "cell_ID")
  
  if(cluster_column!="celltype" & "celltype" %in% names(met)){
    met[["celltype"]] <- NULL
  } 
  data.table::setnames(met, cluster_column, "celltype")
  
  fclist_thresh <-  
  lapply(thresholds, function(thresh){
     incl_cells_thresh <- met[(best_score > thresh & celltype %in% clusters) | is.na(best_score)][["cell_ID"]] 
     clusterwise_foldchange_metrics(normed = normed[markers,incl_cells_thresh]
                                    ,metadata = met[cell_ID %in% incl_cells_thresh]
                                    ,clustercol = "celltype"
                                    )[,thresh:=thresh]
  })
  fcthresh <- rbindlist(fclist_thresh)
  
  nthresh <- fcthresh[,uniqueN(thresh)]
  threshcls <- (colorRampPalette(RColorBrewer::brewer.pal(9, "Reds"))(nthresh))
  names(threshcls) <- fcthresh[,sort(unique(thresh))]
  
  pd <- fcthresh[gene %in% markers][cluster %in% clusters]
  pd[,thresh:=factor(thresh, levels=sort(unique(thresh)))]
  if(!is.null(cluster_order)){
    pd[,cluster:=factor(cluster, levels = cluster_order)]
  }
  if(!is.null(marker_order)){
    pd[,gene:=factor(gene, levels = marker_order)]
  }
  
  p <-  
  ggplot(pd
         ,aes(y = gene, fill=thresh, x = cluster_prop)) + 
    theme_bw() + 
    scale_fill_manual(values = threshcls, guide=guide_legend(reverse=TRUE)) + 
    scale_x_continuous(n.breaks = 14) + 
    geom_bar(stat='identity', position=position_dodge2(),color='black',lwd=0.1) + 
    facet_wrap(~cluster)
  return(p)
}


