#' @title pairsplot for metagene scores and classifications using GGally and ggplot2
#' @description
#' Make a pairsplot for metagene scores and classifications using GGally and ggplot2.
#'  
#' @param metagenes a metagenes object created using `fit_metagene_scores()`
#' @param var one of 'ypost', 'yhat', or 'y' to be plotted (reflecting posterior-averaged, model-predicted, or raw scores of the index gene, respectively).
#' @param nsamples if NULL, plots all cells.  Otherwise, plots a random sample of `nsamples` cells (quicker for large datasets).
#' @param wts optional cells-length vector of wts, reflecting prior probability that cell belongs to **any** of the possible metagene clusters.
#' @param wtlim optional lower weight limit for plotting cells.  If `wts` are provided, limit the plot to cells where `wts > wtlim`.
#' @param colorvar  Optional categorical variable by which to color the cells.  This should be provided as a .
#' @param ptalpha  alpha value controlling transparency of plotted points.
#' @param ptsize  controls size of plotted points.
#' @param obs optional cell-level data frame / data.table that is inner-joined onto the metagene scores.  
#'        This can be used to subset the cells plotted (by only passing `obs` for the cells you want to plot) and/or including `colorvar` variable to be used for plotting.

#' @export
metagene_pairsplot <- function(metagenes, var, nsamples = NULL, wts = NULL
                               ,colorvar=NULL, wtlim = 0.5
                               ,ptalpha = 0.5, ptsize = 0.1, obs = NULL){
  pd <- 
    data.table::rbindlist(
      lapply(names(metagenes), function(xx){
        data.table::data.table(yhat=metagenes[[xx]]$yhat, y=metagenes[[xx]]$y
                   , ypost=metagenes[[xx]]$ypost
                   ,cell_ID = names(metagenes[[xx]]$yhat))[,ct:=xx]
      })
      
    )
  pd <- data.table::dcast(pd, cell_ID ~ ct, value.var=var)
  if(!is.null(wts)){
    pd <- pd[wts > wtlim]
  }
  if(!is.null(obs)){
    pd <- merge(pd, obs, by="cell_ID")
  }
  if(is.null(nsamples)) nsamples <- nrow(pd)
  if(is.null(colorvar)){
    p <- GGally::ggpairs(pd[sample(1:.N, nsamples)][,-c("cell_ID"), with=FALSE]
                         ,lower = list(continuous = GGally::wrap("points", alpha = ptalpha,   size=ptsize))
                         ) + 
      ggplot2::labs(title=var)
      
  } else {
    p <- GGally::ggpairs(pd[sample(1:.N, nsamples)][,-c("cell_ID"), with=FALSE]
                         ,mapping = ggplot2::aes(color=.data[[colorvar]])
                         ,lower = list(continuous = GGally::wrap("points", alpha = ptalpha,   size=ptsize))
                         ) + 
      labs(title=var)
    
  }
  return(p)
}
