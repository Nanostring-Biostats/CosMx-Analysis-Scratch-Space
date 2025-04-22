



#' Plot celltypes in x/y space
#' 
#' @param cluster_column  
#' @param x_column name of column in metadata for x-coordinate
#' @param y_column name of column in metadata for y-coordinate
#' @param cls optional (possibly named) vector of colors to be used
#' @param clusters which clusters to plot.  If NULL (default), then all clusters will be plotted.
#' @param metadata data.frame or data.table of metadata containing cluster_column, x_column, y_column
#' @param ptsize point size for plotting
#' @param alphasize alpha (transparency) for plotting
#' @param plotfirst optional vector of celltypes which should be plotted first (under) the rest.
#' 
#' 
#'
#' @export
#' 
xyplot <- function(cluster_column, x_column = "x_slide_mm", y_column = "y_slide_mm", cls=NULL
                   ,clusters =NULL
                   ,metadata,ptsize=0.25
                   ,plotfirst = NULL
                   ,alphasize=1){
  pd <- data.table::copy(data.table::data.table(metadata))
  if(is.null(clusters)) clusters <- unique(pd[[cluster_column]])
  if(!is.null(plotfirst)){
    plotfirst <- intersect(clusters, plotfirst)
    notplotfirst <- setdiff(clusters, plotfirst)
    p <-ggplot(pd[pd[[cluster_column]] %in% plotfirst]
               ,aes(.data[[x_column]], .data[[y_column]], color=.data[[cluster_column]])) + 
      geom_point(size=ptsize)
    p <- p + geom_point(data=pd[pd[[cluster_column]] %in% notplotfirst]
                        ,aes(.data[[x_column]], .data[[y_column]], color=.data[[cluster_column]])
                        ,size=ptsize, alpha=alphasize) + 
      theme_bw()
    
  } else {
    p <- 
      ggplot(pd[pd[[cluster_column]] %in% clusters]
             ,aes(.data[[x_column]], .data[[y_column]], color=.data[[cluster_column]])) + 
      geom_point(size=ptsize, alpha=alphasize) + 
      theme_bw()
    
  }
  if(is.null(cls)){
    cls <- rep(unname(pals::alphabet()), 100)
    if(any(is.na(suppressWarnings(as.numeric(as.character(clusters)))))){
      clnames <- sort(clusters)
    } else {
      clnames <- sort(as.numeric(as.character(clusters)))
    }
    cls <- cls[1:length(clusters)]
    names(cls) <- clnames
    p <- p + 
      scale_color_manual(values=cls
                         ,guide=guide_legend(override.aes=list(size=4,alpha=1))) #+ 
    
  } else {
    p <- p + 
      scale_color_manual(values=cls
                         ,guide=guide_legend(override.aes=list(size=4,alpha=1))) #+ 
    
  }
  return(p) 
}

