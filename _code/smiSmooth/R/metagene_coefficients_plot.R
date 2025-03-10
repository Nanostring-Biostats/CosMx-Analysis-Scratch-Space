
#' Plot the regression coefficients of metagene scores
#' @param metagenes metagenes object created by fit_metagene_scores
#' @param label_cutoff_abscoef minimum absolute value of regression coefficient to be labeled with geom_label_repel
#' @param label_cutoff_n maximum number of predictors to annotate with geom_label_repel
#' @param max.overlaps max.overlaps argument passed to ggrepel::geom_label_repel for labeling predictor names .
#' @param annotate_nonzero if TRUE, create a text annotation showing the fraction of non-zero predictors out of all possible.
#' 
#' @export
#' 
metagene_coefficients_plot <- function(metagenes
                                       ,label_cutoff_abscoef = 0.05
                                       ,label_cutoff_n = 100
                                       ,max.overlaps = 100
                                       ,annotate_nonzero = TRUE
                                       ){
  
  bhats <- 
    rbindlist(
      lapply(names(metagenes), function(xx){
        as.data.table(metagenes[[xx]]$bhat, keep.rownames=TRUE)[,metagene:=xx][
          order(metagene,-V1),rnk:=1:.N,by=metagene][,.(g=rn,bhat=V1,metagene,rnk)]
      })
    )[order(metagene,rnk)] 

 
 
  plist <-  
    lapply(split(bhats, by = "metagene"), function(xx){
     
      lbld <- xx[abs(bhat) > label_cutoff_abscoef]
      if(nrow(lbld) > label_cutoff_n){
        lbld <- lbld[order(-abs(bhat))][1:label_cutoff_n]
      }
      nz <- xx[,sum(bhat!=0)]
      p <-  
      ggplot(xx, aes(rnk, bhat, label=g)) + geom_point() + 
        facet_wrap(~metagene) + 
        geom_label_repel(data=lbld,max.overlaps = max.overlaps) + 
        geom_hline(yintercept = 0,lty=2,color='red') + 
        theme_bw() + 
        ylab(expression(hat(beta))) + 
        xlab("predictor rank") + 
        annotate("text", x = Inf, y = Inf
                 ,label = paste0(nz, "/", nrow(xx), " non-zero predictors")
                 ,hjust=1.2,vjust = 1.5)
        facet_wrap(~metagene, scales='free') 
      return(p)
    })
  
  return(plotlist = plist)  
}

