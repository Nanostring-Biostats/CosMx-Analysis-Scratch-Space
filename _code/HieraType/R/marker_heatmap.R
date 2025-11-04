
#' Make a marker gene heatmap showing fold changes and proportion of cells expressing the marker gene by cluster.
#' 
#' @param foldchange_metrics metrics calculated using `clusterwise_foldchange_metrics()` function.
#' @param extras if provided, these extra genes will be included in the plot, even if they are not in the `topn` in terms of 
#' `foldchange` for a particular cluster based on fold change or don't meet the `cluster_prop_min` requirement.
#' @param topn plot the `topn` genes (in terms of fold_change, which are expressed in at least cluster_prop_min pct of cells) from each cluster.
#' @param cluster_prop_min minimum proportion of cells that the gene must be expressed in, in order to qualify to be plotted.
#' @param fold_change_min minimum fold change, order to qualify to be plotted.
#' @param featsuse if provided, these are the (only) genes that will be used , in addition to any `extras`
#' @param geneorder if provided, this allows to manually control the order of genes.  
#' Probably should only used if the exact set of features is also provided with `featsuse`.
#' @param clusterorder if provided this allows to manually control the order in which clusters are plotted.
#' @param colorordervar column in foldchange_metrics used in ComplexHeatmap::Heatmap function, which
#'                      determines row and column positions for clusters and genes. 
#' @param orient_diagonal Orient marker genes from bottom-left to top-right (usually easier to read).
#' @param topleft_to_bottomright if TRUE, and orient_diagonal also TRUE, do TL to BR instead of BL to TR.
#'                      
#' @export
#'                      
marker_heatmap <- function(foldchange_metrics
                      ,extras = c()
                      ,featsuse =NULL
                      ,topn=5
                      ,geneorder = NULL
                      ,clusterorder = NULL
                      ,fold_change_min = -1
                      ,cluster_prop_min = 0.05
                      ,colordervar = "scaled_fold_change"
                      ,orient_diagonal = TRUE
                      ,topleft_to_bottomright = FALSE
){
  fctable <- data.table::copy(foldchange_metrics)
  fctable[,scaled_fold_change:=fold_change/max(fold_change),by=gene]
  fctable[,scaled_expr:=cluster_expr/max(cluster_expr),by=gene]
  fctable[is.na(scaled_fold_change),scaled_fold_change:=0]
  fctable[is.na(scaled_expr),scaled_expr:=0]
  
  dmat <- data.table::dcast(fctable, gene ~ cluster, value.var=colordervar)
  fc_mat <- as.matrix(dmat[,-c("gene"),with=FALSE])
  rownames(fc_mat) <- dmat[,gene]
  if(missing(featsuse)){
    featsuse <- unique(c(fctable[order(-fold_change)][fold_change > fold_change_min][cluster_prop > cluster_prop_min,head(.SD,topn),by=cluster][,unique(gene)]
                         , extras))
    
  }
  
  suppressMessages(
    hmap <- 
    ComplexHeatmap::Heatmap(fc_mat[featsuse,])
    )
  suppressMessages(
    grid::grid.grabExpr(hmap <- ComplexHeatmap::draw(hmap))
    )
  if(missing(geneorder)){
    gene_order <- featsuse[ComplexHeatmap::row_order(hmap)]
  } else {
    gene_order <- geneorder
  }
  if(missing(clusterorder)){
    ct_order <- colnames(fc_mat)[ComplexHeatmap::column_order(hmap)]
  } else {
    ct_order <- clusterorder
  }
  pd <- fctable[gene %in% featsuse]
  pd[,gene:=factor(gene, levels = gene_order)]
  pd[,cluster:=factor(cluster, levels = ct_order)]
 
  if(topleft_to_bottomright){
    pd[,cluster:=factor(cluster, levels = rev(ct_order))]
  } 
  
  if(orient_diagonal){
    pd[,rnk:=rank(-fold_change),by=gene]  
    pd[,gene:=factor(gene, levels=pd[rnk==1][order(cluster,-fold_change)][,unique(gene)])] 
  }
  
  hmp <-   
  ggplot(pd, aes(y= cluster, x = gene
                 ,fill=scaled_expr
                 ,size = cluster_prop
  )
  ) + 
    geom_point(pch=21) + 
    scale_fill_gradientn(colors = RColorBrewer::brewer.pal(9,"Reds"),name = "Mean expression in group") + 
    scale_size_continuous(name = "Fraction of cells in group") + 
    theme_classic() +
    theme(axis.text.x=element_text(angle=90,hjust=1,vjust=0.5,face ="bold")
          ,axis.title.x=element_blank()
    ) + 
    theme(axis.text.y=element_text(face ="bold")
          ,axis.title.y = element_blank()
    ) 
   return(hmp) 
}



