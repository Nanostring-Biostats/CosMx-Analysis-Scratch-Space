make.volcano <- function(dat
                         ,compv
                         ,bonf=NULL
                         ,lblsize=4
                         ,titlval=NULL
                         ,sig_color_left = "#1F78B4" ## nice blue
                         ,sig_color_right = "#E31A1C"## nice red
                         ,fold_change_cutoff=1.5
                         ,text_size=20
                         ,custom_label_targets = c()
                         ,interactive = TRUE
                         ,side_label_bool = TRUE
                         ,limit_labels = F
){
  dat <- dat[!is.na(p.value)]
  stopifnot(any(dat[,p.value < 0.9]))
  data <- data.table::copy(dat)[p.value < 0.9]
  data[p.value==0,p.value:=.Machine$double.xmin]
  data <- data[fold_change > 0]
  data <- data[!is.infinite(fold_change)]
  data[,neglog10_pval:=-log10(p.value)]
  data[,log2_FoldChange:=log2(fold_change)]
  if(is.null(bonf)){
    bonf<- data[,0.05/(uniqueN(target) * uniqueN(contrast))]
  } 
  if(is.null(titlval)){
    titlval <- data[1,contrast] 
  }
  
  data[,`:=`(cl="#80808033",bigeff="no")] 
  data[log2_FoldChange < -log2(fold_change_cutoff) & 
         neglog10_pval > -log10(bonf),`:=`(
           cl=sig_color_left
           ,bigeff="neg")] 
  data[log2_FoldChange > log2(fold_change_cutoff) & 
         neglog10_pval > -log10(bonf),`:=`(
           cl=sig_color_right
           ,bigeff="pos")] 
  
  data[bigeff=="no" & (neglog10_pval > -log10(bonf)), `:=`(cl = "grey40",bigeff="sig")] 
  data[is.infinite(neglog10_pval),neglog10_pval:= -log10(.Machine$double.xmin)]
  
  
  # Plot
  max.y <- nDigits(min(data[,min(p.value)], signif(bonf / 2, 1)))
  
  ybrk <- 1 * 10^pretty( c(-1:(max.y)) )
  ylab <- ybrk
  
  clval <-data[,.(cl, bigeff)][,unique(.SD)][,cl]
  names(clval) <-data[,.(cl, bigeff)][,unique(.SD)][,bigeff]
  
  labels <- data[,target]#rownames(data)
  xrng <- data[!is.na(log2_FoldChange)
               ,max(log2_FoldChange) - min(log2_FoldChange)] 
  xrng <- data[!is.na(log2_FoldChange),c(min(log2_FoldChange),max(log2_FoldChange))]
  xrng <- pretty(xrng)
  
  subset_data_for_labeling <- data[p.value < bonf*0.9 | target %in% custom_label_targets]
  subset_data_for_labeling$label <- subset_data_for_labeling[,target]#rownames(subset_data_for_labeling)
  subset_data_for_labeling[fold_change > 1,ndg:= max(xrng) - log2_FoldChange]
  subset_data_for_labeling[fold_change < 1,ndg:= min(xrng) - log2_FoldChange]
  if(limit_labels){
    if(subset_data_for_labeling[,.N] > 50){
      subset_data_for_labeling <- subset_data_for_labeling[abs(log2_FoldChange) > log2(fold_change_cutoff) | target %in% custom_label_targets]
    }
  }
  
  if(interactive == TRUE){
    data$CI <- apply(data, 1, function(x){paste(log2_ci(fold_change_est = as.numeric(x[["fold_change"]]), 
                                                  fold_change_se = as.numeric(x[["SE"]])), collapse = ":")})
    
    vplot <- plotly::plot_ly(type = "scatter", mode = "markers",
                     data, x = ~log2_FoldChange, y = ~p.value,
                     # Hover text:
                     text = ~paste("Target:", target,
                                   "<br>log2 Fold Change: ", round(log2_FoldChange,3),
                                   '<br>95% CI:', CI,
                                   '<br>pval:', prettyNum(p.value, format = "e", digits = 3)),
                     color = ~bigeff,
                     colors = clval) %>%
      plotly::layout( title = titlval,
                      font = list(size = text_size),
                      showlegend = FALSE,
                      yaxis = list(title = 'pvalue', autorange="reversed", 
                                   type = "log", tickformat = '.0e',
                                   tickvals = ybrk, ticktext = formatC(ylab, format = "e", digits = 0)),
                      #xaxis = list(title = paste0(compv[2], " <-- log2(Fold Change) --> ", compv[1])),
                      xaxis = list(title = paste0(compv[2], " "
                                                  , sprintf('\u2190')
                                                  , " log", sprintf('\u2082')
                                                  , "(FC) ", sprintf('\u2192'), " ", compv[1])),
                      margin = list(l = 25, r = 25,
                                    b = 100, t = 75)) %>%
      plotly::add_segments(x = -log2(fold_change_cutoff), xend = -log2(fold_change_cutoff), 
                           y = 0, yend = max(data$p.value)*1.05,
                           line = list(dash = "dash", color = 'black'),inherit = FALSE, showlegend = FALSE) %>%
      plotly::add_segments(x = log2(fold_change_cutoff), xend = log2(fold_change_cutoff), 
                           y = 0, yend = max(data$p.value)*1.05,
                           line = list(dash = "dash", color = 'black'),inherit = FALSE, showlegend = FALSE) %>%
      plotly::add_segments(x = min(data$log2_FoldChange)*1.05, xend = max(data$log2_FoldChange)*1.05, 
                           y = bonf, yend = bonf,
                           line = list(dash = "dash", color = 'black'),inherit = FALSE, showlegend = FALSE)
  }else{
    vplot <- ggplot2::ggplot(data
                             , ggplot2::aes(x=log2_FoldChange
                                            ,y=p.value
                                            ,color=bigeff))  + 
      ggplot2::geom_point() + 
      ggplot2::theme_bw(base_size = text_size)  + 
      ggplot2::scale_color_manual(values=clval,labels=names(clval)) + 
      ggrepel::geom_text_repel(data = subset_data_for_labeling
                               ,ggplot2::aes(log2_FoldChange
                                             ,p.value
                                             ,label = label)
                               ,nudge_x = if(side_label_bool) subset_data_for_labeling[, ndg] else 0
                               ,direction = if(side_label_bool) "y" else c("both", "y", "x")
                               ,max.overlaps = if(limit_labels) 100 else 1000 
                               ,inherit.aes = FALSE
                               ,size=lblsize) + 
      ggplot2::theme(legend.position="none") + 
      ggplot2::geom_hline(yintercept = bonf,linetype=2,color='black') + 
      ggplot2::geom_vline(xintercept = -log2(fold_change_cutoff),linetype=2,color='black') + 
      ggplot2::geom_vline(xintercept = log2(fold_change_cutoff),linetype=2,color='black') + 
      ggplot2::labs(x = paste0(compv[2], " ", sprintf('\u2190'), " log", sprintf('\u2082'), "(FC) ", sprintf('\u2192'), " ", compv[1]) 
                    ,y="pvalue"
                    ,title=titlval) + 
      ggplot2::scale_y_continuous(trans=change_axis_revlog_trans(base=10),breaks=ybrk, 
                                  labels = formatC(ylab, format = "e", digits = 0)) +
      ggplot2::scale_x_continuous(breaks = xrng, limits = range(xrng)) 
  }
  
  return(vplot) 
}

#' change_axis_revlog_trans
#'
#' reverse log transform axis; used to return pvalue rather than -log10(pvalue) on yaxis
#' @param base base in which logs are computed
#' @return revlog_trans reverse log transformation
#' @noRd
change_axis_revlog_trans <- function(base=exp(1)){
  trans <- function(x) -log(x, base)
  inv <- function(x) base^(-x)
  revlog_trans <- trans_new(name=paste0("revlog-", base), 
                            transform=trans, 
                            inverse=inv, 
                            breaks=log_breaks(base=base),
                            domain=c(1e-310, Inf))
  
  return(revlog_trans)
}


#' Make a volcano plot from smide_results
#' 
#' @param smide_results object of class "smide" as returned from `smi_de` or `run_de` functions. 
#' @param comparison one of c("pairwise", "one.vs.rest", "one.vs.all"). 
#'                    See ?smiDE::results for details on the individual types of output tables.
#' @param variable fixed effect variable in DE model for which to make volcano plot.
#' @param targets targets to plot.  typically should be left as the default "all" for full volcano plot. 
#' @param bonf  Optional multiple testing threshold for horizontal dashed line, set by default to 0.05/(# of genes)
#' @param lblsize Size of gene labels for significant genes, passed to geom_text_repel.
#' @param titlval Optional character string for title.
#' @param sig_color_left Color for significant genes with positive negative log fold change.
#' @param sig_color_right Color for significant genes with positive log fold change.
#' @param fold_change_cutoff fold change cutoff for vertical dashed line in plot indicating a 'large magnitude' of effect, used for indicating whether to color targets as non-grey.
#' @param text_size size of plot text set using ggplot2::theme(text=element_text(size=text_size))
#' @param interactive TRUE = plotly, FALSE = static ggplot
#' @param side_label_bool TRUE = labels nudged to the side of the plot.  FALSE = labels places at datapoint (with lines as determined by ggrepel)
#' @param limit_labels FALSE = For non-interactive plots, plot labels for all significant genes with p-value < 0.9*bonf threshold (ggrepel::geom_text_repel max.overlaps argument is set to 1,000).  
#'                     TRUE = genes with fold_changes with lesser magnitude than `fold_change_cutoff` will not be labeled, and ggrepel::geom_text_repel max.overlaps argument will be limited to 100 genes.
#'                     
#' @return A ggplot2 (volcano) plot object, with log2 fold changes plotted on x-axis, and 
#'
#' @examples
#' library(cowplot)
#' datadir <- system.file("extdata", package="smiDE")
#' de_results <- readRDS(paste0(datadir, "/sample_de_results.rds"))
#' 
#' dat <- results(de_results, comparisons="one.vs.rest", variable="group")[["one.vs.rest"]]
#' head(dat)
#' volcano(de_results, "one.vs.rest", "group", interactive = TRUE)
#' volcano(de_results, "pairwise", "group", interactive = FALSE)
#' plist <- volcano(de_results, "pairwise", "group",text_size=8, lblsize=2, interactive = FALSE)
#' names(plist)
#' cowplot::plot_grid(plotlist = plist)
#' 
#' @export
volcano <- function(smide_results
                    ,comparison = "one.vs.rest"
                    ,variable = NULL
                    ,targets = "all"
                    ,bonf=NULL
                    ,lblsize=4
                    ,titlval=NULL
                    ,sig_color_left = "blue"
                    ,sig_color_right = "red"
                    ,fold_change_cutoff=1.5
                    ,text_size=20
                    ,custom_label_targets = c()
                    ,interactive = TRUE
                    ,side_label_bool = FALSE 
                    ,limit_labels = FALSE
){
  
 
  if(inherits(smide_results, "smide")){
    if(missing(variable)) variable <- smide_results$groupVar
    if(is.null(comparison)) stop("comparison argument (one of 'pairwise', 'one.vs.rest', or 'one.vs.all')\nmust be specified when passing in smiDE results object.")
    dat <- results(smide_results
                   ,comparisons = comparison[1]
                   ,variable
                   ,targets
    )[[comparison[1]]]
  } else {
    dat <- data.table::copy(smide_results)
    if(!all(c("contrast", "p.value", "fold_change", "SE", "target") %in% names(dat))){
      misscols <- paste0(setdiff(c("contrast", "p.value", "fold_change", "SE", "target"), names(dat))
             ,collapse=", ")
      
      stop(paste0(" required columns: ", misscols, " not found.\nInput datatable typically created using `smiDE::results()` function."))
    }
  }
  
  
  plotlist <- list()
  for(cc in dat[!is.na(contrast),unique(contrast)]){
    
    pdat <- dat[contrast==cc]
    compv <- strsplit(as.character(cc), "\\ vs.\\ ")[[1]]
    if(grepl(" [\\/-]\\ ", cc)) compv <- strsplit(as.character(cc), "\\ [\\/-]\\ ")[[1]]
    
    plotlist[[cc]] <- make.volcano(pdat
                                   ,compv
                                   ,bonf
                                   ,lblsize
                                   ,titlval
                                   ,sig_color_left
                                   ,sig_color_right
                                   ,fold_change_cutoff
                                   ,text_size
                                   ,custom_label_targets
                                   ,interactive
                                   ,side_label_bool = side_label_bool
                                   ,limit_labels = limit_labels
    )
    
  }
  return(plotlist)
}

#' calculate CI of fold chage
#'
#' @param fold_change_est estimated fold change value 
#' @param fold_change_se standard error of fold change
#' @return 95% confidence interval of fold change 
#' @noRd
log2_ci <- function(fold_change_est, fold_change_se){
  l2_est <- log2(fold_change_est)
  ### Delta method
  l2_se <- sqrt((1/(fold_change_est*log(2)))^2 * fold_change_se^2)
  return(round(l2_est + c(-1, 1)*qnorm(0.975)*l2_se, 4))
}


#' Save volcano plot files from list
#'
#' @param volcanoPlot list of volcano plot figures
#' @param outFolder output directory 
#' @examples 
#' \dontrun{
#' 
#' datadir <- system.file("extdata", package="smiDE")
#' de_results <- readRDS(paste0(datadir, "/sample_de_results.rds"))
#' 
#' dat <- results(de_results, comparisons="one.vs.rest", variable="group")[["one.vs.rest"]]
#' head(dat)
#' 
#' plist <- volcano(de_results, "one.vs.rest", "group", interactive = TRUE)
#' # tempDir <- tempdir()
#' # saveVolcano(volcanoPlot = plist, outFolder = getwd())
#' 
#' }
#' 
#' @export
saveVolcano <- function(volcanoPlot, outFolder){
  for(i in names(volcanoPlot)){
    name <- gsub("\\W", "-", i)
    name <- gsub("_", "-", name)
    if(class(volcanoPlot[[i]])[2] == "ggplot"){
      ggplot2::ggsave(filename = paste0(outFolder, "/", name, ".png"), 
                      plot = volcanoPlot[[i]], device = "png", width = 7, height = 7)
    } else if(class(volcanoPlot[[i]])[2] == "htmlwidget"){
      htmlwidgets::saveWidget(
        widget = volcanoPlot[[i]], # the plotly object
        file = paste0(outFolder, "/", name, ".html")
      )
      unlink(paste0(outFolder, "/", name, "_files"), recursive = TRUE)
    } else {
      stop("Invalid volcano plot given")
    }
  }
}

#' return number of digits in a number
#'
#' @param x number 
#' @return number of digits in x
#' @noRd
nDigits <- function(x) suppressWarnings(as.integer(strsplit(x = formatC(formatC(x = x, format="e"), format = "s"), split = "e")[[1]][2]))

