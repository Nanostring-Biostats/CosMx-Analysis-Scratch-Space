

## Code for plots from NSCLC DE Analysis:
## Fig 3 (Macrophage cell DE: Segmentation error example, DAG, fold change with/without attenuation, volcano plot)
## Supp Fig 5  Significant genes and annotations from Human Protein Atlas by adjustment strategy


cast_and_melt <- function(dat, castvar, valuevar, idvars, reflevel){
  dw <- dcast(dat, as.formula(paste0(paste0(idvars, collapse = "+"), "~", paste0(castvar,collapse="+"))), value.var=valuevar)
  colname_w <- colnames(dw)
  dm <- 
    tryCatch({
      melt(dw, id.vars=c(idvars, reflevel))
    }, error = function(ee){
      message(paste0(colname_w, collapse=",")) 
    })
  return(dm)              
}

library(data.table); setDTthreads(2)
library(ggtext)
library(ggrepel)
library(RColorBrewer)

ms <- fread("nsclc_de_results/ms.plotting.csv")

ovrplotuse <- fread("nsclc_de_results/ovr.plotting.csv")
ovrplotuse <- ovrplotuse[term=="niche"][contrast!=""][
  ,c("contrast", "cell_type_label", "frmla", "slide_use"
     ,"groupadj","fam"
     ,"fold_change", "p.value"
     , "imm_anno"
     ,"maxratio","cell_type_n", "cell_type_cts","cpc_cell_type", "ncells_1", "counts_1", "cpc_1", "target")
  ,with=FALSE
]
ovr <- copy(ovrplotuse)
vd <- ovr[frmla %in% c("naiive", "contam_otherct")][
  slide_use=="all"][
    groupadj=="ranef_tissue"][
      fam=="nbinom2"]

vd[,hpar_label:=gsub("^.*\\ enriched","enhanced/enriched", imm_anno)]
vd[,hpar_label:=gsub("^.*\\ enhanced","enhanced/enriched", hpar_label)]
vd[grepl("monocyte|macrophage|myeloid|DC",hpar_label), hpar_label:=gsub("\\(.*\\)", "(monocyte,macrophage,myeloid,DC)",hpar_label)]
vd[grepl("T-cell|Treg |T-reg",hpar_label), hpar_label:=gsub("\\(.*\\)", "(T cell)",hpar_label)]
vd[grepl("plasma|B",hpar_label), hpar_label:=gsub("\\(.*\\)", "(B cell)",hpar_label)]
vd[,hpar_match:=hpar_label]
vd[grepl("^T",cell_type_label) & grepl("T cell", hpar_label),hpar_match:="enhanced_in_cell_type"]
vd[grepl("macrophage",cell_type_label) & grepl("macro", hpar_label),hpar_match:="enhanced_in_cell_type"]
vd[!grepl("Low|Not|enhanced_in_", hpar_match),hpar_match:="enhanced/enriched (another immune cell type)"]
vd[is.na(hpar_label), hpar_match:="no annotation"]

nsig <- 
  vd[cell_type_label %in% c("T CD4", "T CD8", "Treg", "macrophage"),.(nsig=sum(p.value < 0.05/(960*9))
                                                                      ,nsig_filter=sum(p.value < 0.05/(960*9) & maxratio < 1))
     ,by=.(hpar_match
           ,cell_type_label
           ,frmla, fam, groupadj)]

nsigw <- 
  dcast(nsig
        ,hpar_match + cell_type_label~ frmla
        , value.var=c("nsig", "nsig_filter")
        ,fill=0)

mnsig <- melt(nsigw, id.vars=c("hpar_match", "cell_type_label"))

mnsigp <- copy(mnsig)
mnsigp[,hpar_match:=factor(hpar_match, levels=mnsigp[order(-value),unique(hpar_match)])]
mnsigp[,adjustment_strategy:=factor(variable
                                    ,levels=c("nsig_naiive"
                                              ,"nsig_contam_otherct"
                                              ,"nsig_filter_naiive"
                                              ,"nsig_filter_contam_otherct")
                                    ,labels=c("No adjustment"
                                              ,"Covariate Adjusted"
                                              ,"Contamination Ratio Filter"
                                              ,"Covariate Adjusted + Contamination Ratio Filter")
)]


mnsigp[,adjustment_strategy:=factor(variable
                                    ,levels=c("nsig_naiive"
                                              ,"nsig_contam_otherct"
                                              ,"nsig_filter_naiive"
                                              ,"nsig_filter_contam_otherct")
                                    ,labels=c("No adjustment"
                                              ,"Covariate Adjusted"
                                              ,"Overlap Ratio Filter"
                                              ,"Covariate Adjusted + Overlap Ratio Filter")
)]

###############################################################
#######           Supplementary Figure 5
###############################################################
barp_adj_strategy <- 
  ggplot(mnsigp, aes(x=hpar_match, y=value, fill=adjustment_strategy)) + 
  geom_bar(stat='identity',position=position_dodge(width=0.9)) + 
  facet_wrap(~cell_type_label,scales='free_y') + 
  theme_bw() + 
  geom_text(data=mnsigp, aes(x=hpar_match, y=value + 8, fill=adjustment_strategy, label=value)
            ,inherit.aes=FALSE, position=position_dodge(width=0.9)) + 
  coord_flip() + 
  theme(text=element_text(face='bold',size=12))   + 
  labs(title="# of Significant DE genes across spatial niches\nby adjustment strategy and cell type"
       ,y="# of significant genes (p < 0.05/(960 genes x 9 contrasts))"
       ,x="Condensed Gene Label Category from Human Protein Atlas") + 
  theme(legend.position="bottom")
barp_adj_strategy
ggsave("plots/barplot_adjustment_strategy.png"
       ,units='px'
       ,width=1862*sqrt(10)
       ,height=867*sqrt(10)
       ,dpi=350
) 

################################################################################
####                                                             ############### 
####   Figure 3 Segmentation errors and impact on DE             ###############
####                                                             ###############
################################################################################
library(magick)
library(cowplot)
bonf_contrast <- 0.05/(960*9)

pb <- ggdraw() +
  draw_image("nsclc_de_results/DAG_segoverlap.PNG")

pa <- ggdraw() +
  draw_image("nsclc_de_results/krt17_macrophage_contamination.PNG")

cowplot::plot_grid(pa, pb, rel_widths = c(1, 0.6), rel_heights = c(1,0.2), labels='auto', label_size=24)

vdplot <- ovr[frmla %in% c("naiive", "contam_otherct")][slide_use=="all"][groupadj=="ranef_tissue"][fam=="nbinom2"]

ctt <- "macrophage"
bonf_contrast <- 0.05/(960*9)
vdplot[,cell_type:=cell_type_label]
for(ctt in ovr[,cell_type]){
  vd <- vdplot[cell_type==ctt]
  vd[,hpar_label:=gsub("^.*\\ enriched","enhanced/enriched", imm_anno)]
  vd[,hpar_label:=gsub("^.*\\ enhanced","enhanced/enriched", hpar_label)]
  vd[grepl("monocyte|macrophage|myeloid",hpar_label), hpar_label:=gsub("\\(.*\\)", "(monocyte,macrophage,myeloid)",hpar_label)]
  vd[grepl("B-cell",hpar_label), hpar_label:=gsub("\\(.*\\)", "(B cell)",hpar_label)]
  vd[grepl("T-cell|Treg |T-reg",hpar_label), hpar_label:=gsub("\\(.*\\)", "(T cell)",hpar_label)]
  vd[is.na(hpar_label), hpar_label:="no annotation"]
  vd[,contam_filter:="Passes Filter"]
  vd[maxratio > 1, contam_filter:="Fails Filter"]
  names_hpar_cols <- vd[p.value < bonf_contrast,.N,by=hpar_label][order(-N)][,unique(hpar_label)]
  names_hpar_cols <- c(names_hpar_cols, "enhanced/enriched (another immune cell type)")
  hpar_cols <- pals::alphabet(n=length(names_hpar_cols)) 
  hpar_cols[length(hpar_cols)] <- "red"
  names(hpar_cols) <- names_hpar_cols
  bonf_contrast_filter <- vd[maxratio <= 1,0.05/(uniqueN(target) * uniqueN(contrast))]
  
  
  
  pdsplit <- split(vd[contrast %in% c("tumor interior vs. avg.rest", "plasmablast-enriched stroma vs. avg.rest", "lymphoid structure vs. avg.rest")]
                   ,by="contrast")
  naivepl <- list()
  adjpl <- list()
  volcpl <- list()
  cc <- 3
  pd <- pdsplit[[names(pdsplit)[cc]]]
  pd[!grepl("Low|Not|monocyte", hpar_label),hpar_label:="enhanced/enriched (another immune cell type)"]
  pd[,hpar_label:=factor(hpar_label, levels=c("enhanced/enriched (monocyte,macrophage,myeloid)"
                                              ,"Low immune cell specificity"
                                              ,"Not detected in immune cells"
                                              ,"enhanced/enriched (another immune cell type)"))]  
  pd[,ratiocut:=cut(maxratio, c(0, 0.5, 0.8,0.9,1,1.1,1.2,1.5,2,max(maxratio)))]
  
  pd[,neglog10_pval:=-log10(p.value)] 
  pd[is.infinite(neglog10_pval),neglog10_pval:= -log10(.Machine$double.xmin)]
  
  
  pd[,adjustment:=frmla] 
  pd[frmla=="naiive", adjustment:="Unadjusted"]
  pd[frmla=="contam_otherct", adjustment:="Covariate-Adjusted"]
  pd[,adjustment:=factor(adjustment, levels=c("Unadjusted", "Covariate-Adjusted"))]
  
  
  # Plot
  brkpt <- 15 
  nsqz <- 10
  brkpt <- 100 
  nsqz <-  3
  brk.y <- ceiling(pd[-log10(p.value) < brkpt,max(neglog10_pval)])
  max.y <- max(pd[,ceiling(max(neglog10_pval))], ceiling(bonf_contrast + 1))
  sqz <- seq(brk.y, max.y, length.out=nsqz) 
  sqzy <- seq(brk.y, brk.y + 2, length.out=nsqz) 
  sqzyf <- function(np){
    pct <- (np - min(sqz))/max(sqz) 
    return(pct * (max(sqzy) - min(sqzy)) + min(sqzy))
  }
  
  pd[neglog10_pval > brk.y,neglog10_pval:=sqzyf(neglog10_pval)] 
  ybrk <- c(1:(brk.y-1))
  ylab <- round(c(1:(brk.y-1)))
  breakit <- pd[,sum(neglog10_pval > brk.y) > 0]
  if(breakit){
    ybrk <- c(ybrk, sqzy)  
    ylab <- round(c(ylab, sqz))
  } 
  
  
  rmp <- grDevices::colorRampPalette(brewer.pal(8,"RdYlBu"))
  ratiocls <- rev(rmp(pd[,uniqueN(ratiocut)]))
  names(ratiocls) <- pd[,levels(ratiocut)]
  
  xlab <- pd[1,contrast]
  xlab <- strsplit(xlab, " vs. ")[[1]]
  xlab <- bquote(.(xlab[2]) %<-% '     log2(fold change)     ' %->% .(xlab[1]))
  
  
  pdlbl <- pd[abs(log2(fold_change)) > log2(1.2)][p.value < bonf_contrast]
  pdlbl[,hpar_label:=factor(hpar_label, levels=c("enhanced/enriched (monocyte,macrophage,myeloid)"
                                                 ,"Low immune cell specificity"
                                                 ,"Not detected in immune cells"
                                                 ,"enhanced/enriched (another immune cell type)"))]  
  hpar_cols_use <- hpar_cols[levels(pdlbl[["hpar_label"]])]
  
  pdlbl[,nudge_x:=2 - log2(fold_change)] 
  pdlbl[log2(fold_change) < 0,nudge_x:=-2 - log2(fold_change)] 
  rplsize <- 2.5
  naivepl[[cc]] <-  
    ggplot(pd[-log10(p.value) < -log10(bonf_contrast)][frmla=="naiive"][counts_1 > 0.005][p.value < 0.9]
           ,aes(x=log2(fold_change),y=-log10(p.value))) + 
    geom_point() + 
    geom_point(data=pd[-log10(p.value) >= -log10(bonf_contrast)][frmla=="naiive"]
               ,aes(x=log2(fold_change),y=-log10(p.value), color=hpar_label
                    #      ,aes(x=log2(fold_change),y=neglog10_pval, color=hpar_label
                    ,shape=contam_filter)) + 
    theme_bw() + 
    labs(x=xlab) + 
    scale_color_manual(name="HPAR annotation", values=hpar_cols_use
                       ,guide=guide_legend(nrow=2,override.aes=list(size=4))) + 
    scale_shape_manual(name="Overlap Ratio Filter", values=c(4,19), guide=guide_legend(override.aes=list(size=4))) + 
    geom_label_repel(data=pdlbl[frmla=="naiive"]
                     , aes(x=log2(fold_change), y=-log10(p.value), label=target)
                     ,size=rplsize
                     #, aes(x=log2(fold_change), y=neglog10_pval, label=target)
    ) + 
    geom_hline(yintercept = -log10(bonf_contrast), linetype=2,color='red') + 
    theme(legend.position="bottom", legend.box = "vertical")   + 
    labs(title=paste0(pd[1,contrast]))  + 
    theme(text=element_text(size=14, face='bold'))#+ 
  
  
  xycomp <- dcast(pd, contam_filter + ratiocut + maxratio + hpar_label + target ~ frmla, value.var=c("fold_change", "p.value"))
  attenuation_by_ratio <- 
    ggplot(xycomp
           ,aes(x=log2(fold_change_naiive),y=log2(fold_change_contam_otherct)
                ,color=ratiocut
                ,shape=contam_filter
           )) + 
    geom_point() + 
    scale_color_manual(name='Overlap Ratio', values=ratiocls, guide=guide_legend(override.aes=list(size=4,alpha=1))) + 
    scale_shape_manual(name='Overlap Ratio\nFilter',values=c(4,19), guide=guide_legend(override.aes=list(size=4))) + 
    geom_abline(slope=1,intercept=0, linetype=2,color='red') + 
    theme_bw() + 
    geom_hline(yintercept=0, linetype=2,color='red') + 
    geom_vline(xintercept=0, linetype=2,color='red') + 
    labs(title = "Attenuation of fold change adjusting for neighbor cell expression\nand overlap ratio as a filtering metric"
         ,x="Unadjusted log2(fold change)", y = "Covariate-Adjusted log2(fold change)") +
    coord_fixed() + 
    theme(text=element_text(size=14, face='bold')) + 
    theme(plot.title=element_text(size=14)) 
  
  attenuation_by_pvalue <- 
    ggplot(xycomp
           ,aes(x=-log10(p.value_naiive),y=-log10(p.value_contam_otherct)
                ,color=ratiocut
                ,shape=contam_filter
           )) + 
    geom_point() + 
    theme_bw() + 
    scale_color_manual(name='Overlap Ratio', values=ratiocls, guide=guide_legend(override.aes=list(size=4))) + 
    scale_shape_manual(name='Overlap Ratio\nFilter', values=c(4,19), guide=guide_legend(override.aes=list(size=4))) + 
    geom_abline(slope=1,intercept=0, linetype=2,color='red') + 
    geom_hline(yintercept=0, linetype=2,color='red') + 
    geom_vline(xintercept=0, linetype=2,color='red') + 
    coord_fixed()
  adjcoef_xy <- 
    merge(ms[grep("otherct_expr",term)][
      ,.(term, est, se, pval, target, groupadj, cell_type, fam, frmla, slide_use)]
      ,pd
      ,by=c("target", "cell_type", "slide_use","frmla","groupadj" , "fam"))
  
  adj_vs_ratio <-  
    ggplot(adjcoef_xy , aes(x=maxratio, y=est, color=hpar_label
                            #                  , shape=contam_filter
    )
    ) + 
    geom_point(alpha=0.5) + 
    theme_bw() + 
    geom_smooth(method="loess",se=FALSE) + 
    labs(y="Adjustment Regression Coefficient", x="log(Contamination Ratio)") + 
    scale_x_continuous(trans='log', labels = scales::number_format(0.01),n.breaks=20) +
    scale_color_manual(name="HPAR annotation", values=hpar_cols_use
                       ,guide=guide_legend(nrow=2,override.aes=list(size=4))) + 
    theme(legend.position = "none")
  
  xlab <- pd[1,contrast]
  xlab <- strsplit(xlab, " vs. ")[[1]]
  xlab <- gsub("^", "  ", xlab) 
  xlab <- gsub("$", "  ", xlab) 
  xlab <- bquote(.(xlab[2]) %<-% '   log2(fold change)' %->% .(xlab[1]))
  
  adjpl[[cc]] <-  
    ggplot(pd[-log10(p.value) < -log10(bonf_contrast)]
           ,aes(x=log2(fold_change),y=-log10(p.value))) + 
    geom_point() + 
    geom_point(data=pd[-log10(p.value) >= -log10(bonf_contrast)]
               ,aes(x=log2(fold_change),y=-log10(p.value), color=hpar_label
                    ,shape=contam_filter)) + 
    theme_bw() + 
    facet_wrap(~adjustment) + 
    labs(x=xlab,title='macrophages: genes DE in tumor-interior') + 
    theme(legend.position="bottom"
          #,legend.box = "vertical"
          ,legend.box.just = "left"
          ,legend.spacing.y=unit(1,'lines')
          ,legend.title=element_text(margin = margin(r=50))
          #,legend.margin=margin(0,0,0,0)
    ) + 
    scale_color_manual(name="HPAR annotation", values=hpar_cols_use
                       ,guide=guide_legend(nrow=4,override.aes=list(size=4)
                       #,guide=guide_legend(nrow=2,override.aes=list(size=4)
                                           ,legend.spacing.y=unit(0, 'lines')
                                           ,legend.margin=margin(0,0,0,0)
                       )) + 
    scale_shape_manual(name="Overlap Ratio Filter", values=c(4,19), guide=guide_legend(override.aes=list(size=4)
                                                                             #   ,nrow=1
                                                                                ,nrow=2
                                                                                ,legend.spacing.y=unit(0, 'lines')
                                                                                ,legend.margin=margin(0,0,0,0)
    )
    ) + 
    geom_label_repel(data=pdlbl, aes(x=log2(fold_change), y=-log10(p.value)
                                     ,label=target
                                     ,color=hpar_label)
                     ,size=rplsize
                     #  , max.overlaps=100
                     , show.legend = FALSE) + 
    geom_hline(yintercept = -log10(bonf_contrast), linetype=2,color='red') + 
    theme(text=element_text(size=14, face='bold'))
  
  
  pbd <- cowplot::plot_grid(pb, adjpl[[cc]]
                            ,rel_heights = c(0.35,0.65), nrow=2, labels=c("b", "d"), label_size = 24) + 
    theme(plot.background = element_rect(fill = 'white',color='white')) + 
    cowplot::panel_border(remove=TRUE) 
  
  pac <- cowplot::plot_grid(pa, attenuation_by_ratio
                            ,rel_heights = c(0.5,0.5), nrow=2, labels=c("a", "c"), label_size = 24) + 
    theme(plot.background = element_rect(fill = 'white',color='white')) + 
    cowplot::panel_border(remove=TRUE) 
  
  pabcd <-  
  cowplot::plot_grid(pac, pbd, rel_widths = c(0.6,1), ncol=2, labels="")  + 
    theme(plot.background = element_rect(fill = 'white',color='white')) + 
    cowplot::panel_border(remove=TRUE) 
  pabcd

  ##### Figure3  
  ggsave("plots/Figure3.png"
         ,plot = pabcd
         ,units='px'
         ,width=1705*sqrt(30)
         ,height=890*sqrt(30)
         ,dpi=400
  ) 
    
  
}


