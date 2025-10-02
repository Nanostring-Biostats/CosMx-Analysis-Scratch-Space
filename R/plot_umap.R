

#' A UMAP plotting utility function
#' 
#' @description
#' A helper function for plotting UMAP using ggplot2 and ggrepel packages.
#'
#' @param umapreduc name of umap dimension reduction in a `semuse` seurat object. 
#' @param clustercol name of column to color by in `semuse@meta.data`. 
#' @param cls optional (possibly named) vector of colors to plot categories of `clustercol` by.
#' @param plotfirst optional vector of categories in `clustercol` to be plotted first 
#' (potentially allowing other less frequent categories to be highlighted by overlaying on top).
#' @param alpha transparency of points to be plotted
#' @param max.overlaps max.overlaps parameter passed to `ggrepel::geom_label_repel`
#' @param lblsize size of labels to be passed to `ggrepel::geom_label_repel`
#'
#' @export  
plot_umap <- function(umapreduc,clustercol, semuse, cls=NULL,plotfirst=NULL,alpha=1,max.overlaps = 100,lblsize = 4
                      ,cellid_colname = "cell_ID"){
  
  metadata <- data.table::copy(semuse@meta.data)
  metadata[["cell_ID"]] <- rownames(metadata)
   
  umapd <-  
    data.table(semuse@reductions[[umapreduc]]@cell.embeddings
               ,keep.rownames = TRUE)
  setnames(umapd, c(names(umapd)[2:3]), c("UMAP_1", "UMAP_2"))
  obsmrk <- merge(data.table::data.table(metadata), umapd
                  ,by.x="cell_ID"
                  , by.y="rn")
  
  obstxt <- obsmrk[,lapply(.SD, median),by=c(clustercol),.SDcols=paste0("UMAP_",1:2)]
  
  clusters <- unique(obsmrk[[clustercol]])
  
  if(!is.null(plotfirst)){
    p <- 
      ggplot(obsmrk[obsmrk[[clustercol]] %in% plotfirst]
             , aes(UMAP_1, UMAP_2, color=.data[[clustercol]])) + 
      geom_point(size=0.2,alpha=alpha) + 
      geom_point(data=obsmrk[!obsmrk[[clustercol]] %in% plotfirst]
                 ,aes(UMAP_1, UMAP_2, color=.data[[clustercol]]), size=0.2,alpha=alpha) + 
      theme_bw() + coord_fixed() + 
      geom_label_repel(data=obstxt, aes(x=UMAP_1, y=UMAP_2, label=.data[[clustercol]]),show.legend=FALSE
                       ,inherit.aes=FALSE,color='black', max.overlaps = max.overlaps, size = lblsize)
    
  } else {
    p <- 
      ggplot(obsmrk, aes(UMAP_1, UMAP_2, color=.data[[clustercol]])) + 
      geom_point(size=0.2,alpha=alpha) + 
      theme_bw() + coord_fixed() + 
      geom_label_repel(data=obstxt, aes(x=UMAP_1, y=UMAP_2, label=.data[[clustercol]]),show.legend=FALSE
                       ,inherit.aes=FALSE,color='black', max.overlaps = max.overlaps, size = lblsize)
    
  }
  if(is.null(cls)){
    cls <- rep(unname(pals::alphabet()), 10)
    if(any(is.na(suppressWarnings(as.numeric(as.character(clusters)))))){
      clnames <- clusters 
    } else {
      clnames <- sort(as.numeric(as.character(clusters)))
    }
    cls <- cls[1:length(clusters)]
    names(cls) <- clnames
    p <- p +  
      scale_color_manual(values=cls # rep(unname(pals::alphabet()), 3)
                         ,guide=guide_legend(override.aes=list(size=4,alpha=1)))
    
  } else {
    p <- p +  
      scale_color_manual(values=cls
                         ,guide=guide_legend(override.aes=list(size=4,alpha=1)))
  }
  return(p)
}