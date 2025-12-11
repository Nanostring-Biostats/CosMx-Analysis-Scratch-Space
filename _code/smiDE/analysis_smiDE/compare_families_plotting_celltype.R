
##
## Code for plots from NSCLC DE Analysis:
## Fig 2 (Macrophage cell DE: effect size, signed -log10 p values, computation time comparison, venn diagram)
## Supp Fig 1 (effect sizes by family)
## Supp Fig 2 (signed -log10 p values by family)
## Supp Fig 3 (venn diagrams of overlapping significant genes) 


library(data.table); setDTthreads(2)
library(ggvenn)
library(ggVennDiagram)
library(ggrepel)
library(RColorBrewer)


ovrplotuse <- fread("nsclc_de_results/ovr.plotting.csv")
ovrplotuse <- ovrplotuse[frmla=="contam_otherct"][term=="niche"][contrast!=""][
  ,c("contrast", "cell_type_label", "frmla", "slide_use"
     ,"groupadj","fam"#, "sand_se"
     ,"fold_change", "p.value"
     , "imm_anno"
     ,"maxratio","cell_type_n", "cell_type_cts","cpc_cell_type", "ncells_1", "counts_1", "cpc_1", "target")
     ,with=FALSE
  ]


### consistent counts_1, cpc_1 across families
###  - these are different calculated differently for 'gaussian' , which uses total count normalized data
### use raw counts instead
raw_counts_niche <- 
ovr[fam=="poisson"][,c("contrast", "cell_type_label", "frmla", "slide_use"
     ,"counts_1", "ncells_1", "cpc_1","target")
     ,with=FALSE][,unique(.SD)]

ovrplotuse <- 
merge(raw_counts_niche, ovrplotuse[,-c("counts_1", "ncells_1", "cpc_1"),with=FALSE]
      ,by=c("contrast", "cell_type_label", "frmla", "slide_use","target"))

ovrplotuse[,variablef_m:=paste0(groupadj, "_FALSE_", fam)]
ovrplotuse[variablef_m=="ranef_tissue_FALSE_poisson",  modeltype:="poisson RE"]
ovrplotuse[variablef_m=="ranef_tissue_FALSE_gaussian", modeltype:="gaussian RE"]
ovrplotuse[variablef_m=="ranef_tissue_FALSE_nbinom2",  modeltype:="negative binomial RE"]

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


#### Set colors for plots
fam_cols <- c("#1B9E77", "#D95F02", "#7570B3")
names(fam_cols) <- c("negative binomial","poisson", "gaussian")

cell_type_cols_all <- c(brewer.pal(8,"Set1"))[1:5]
names(cell_type_cols_all) <- c("macrophage", "T CD4", "T CD8","Treg", "tumor")

ovrplotuse <- ovrplotuse[groupadj=="ranef_tissue"]
ctt <- "plasmablast"
palist <- pblist <- pclist <- list()

rplsize <-3 
pdl <- list()
annotl <- list()
for(ctt in unique(ovrplotuse[["cell_type_label"]])){
  ovrplot <- ovrplotuse[cell_type_label==ctt]

  dm_fc <- cast_and_melt(
    ovrplot
    ,castvar = c("fam")
    ,valuevar="fold_change"
    ,idvars = c("contrast", "cell_type_label", "frmla", "groupadj", "slide_use"
                ,"target", "imm_anno", "maxratio","cell_type_n"
                ,"ncells_1", "counts_1", "cpc_1", "cpc_cell_type", "cell_type_cts")
    ,reflevel = "nbinom2"
  )
  
  pd <- dm_fc[!is.na(value)][!is.na(nbinom2)]
  pd[value < 0,value:=0]
  pd[,lbl:=paste0(cell_type_label, " (", scales::comma(cell_type_n), " cells)")]
  pd[,lbl:=factor(lbl, levels=pd[order(cell_type_n)][,unique(lbl)])]
  pd[,ncut:=cut(cell_type_n, c(0, 1000, 5000, 10000, 50000, 10**6)
                ,labels=c("<1k cells", "1-5k cells", "5-10k cells", "10-50k cells", ">50k cells"), include.lowest=TRUE)]
  pd[,.N,by=.(cell_type_n, lbl,ncut)]
  pd[,ncut:=cut(cell_type_n, c(0, 5000, 10000, 50000, 200000, 10**6)
                ,labels=c("<5k cells", "5-10k cells", "10-50k cells", "50-200k cells", "227k cells"), include.lowest=TRUE)]
  
  pd[,summary(cpc_cell_type)]
  pd[,cpc_cut:=cut(cpc_cell_type, c(0, 0.01, 0.05, 0.1, 100)
                   ,labels = c("< 0.01", "0.01-0.05", "0.05-0.1", "> 0.1")
                   , include.lowest=TRUE)]
  
  
  #pd[maxratio < 1][cpc_cell_type > 0.01][cts > 10][log2(value) < -4][cell_type_label=="tumor"]
  pd[maxratio < 1][cpc_cell_type > 0.01][counts_1 > 10][variable=="gaussian"][cell_type_label=="tumor"][value ==0]
  annot <- pd[maxratio < 1][cpc_cell_type > 0.01][counts_1 > 10][value!=0,.(annotx=quantile(log2(nbinom2),0.5)
                                                                    ,annoty=quantile(log2(value), 0.0001)
                                                                    #,bhat=coef(lm(log2(value) ~ log2(nbinom2)-1)[1])
                                                                    ,rhohat=cor(log2(value),log2(nbinom2))
                                                                    #,
                                                                    )
                                                                 ,by=.(variable,lbl, cell_type_label)]
  
  annot[,annotx:=-Inf]
  annot[,annoty:=Inf]
  annot[,hjustvar:=-0.5]
  annot[,vjustvar:=(1:.N)*1.5]

  annotl[[ctt]] <- annot 
  pdl[[ctt]] <- pd 
  pa <- 
  ggplot(pd[maxratio < 1][cpc_cell_type > 0.01][counts_1 > 10][value!=0]#[cpc_cell_type > 0.01]#[cts > 10]#[groupadj=="fixef_tissue"]
         ,aes(x=log2(nbinom2)
              ,y = log2(value)
              ,color=variable
              ,group=variable
              ,shape=cpc_cut
              # ,shape=contrast
         )) + 
    geom_point(size=1) + 
    theme_bw() + 
    scale_color_manual(values=fam_cols, guide=guide_legend(override.aes=list(size=4))
                       ,name="Alternative Distribution") + 
    scale_shape_manual(values=c(4,2,1, 7)
                       ,guide=guide_legend(override.aes=list(size=4))
                       ,name="counts per cell") + 
    geom_abline(slope=1,intercept=0) + 
    geom_smooth(method='lm',se=FALSE)  + 
    labs(x = "log2(fold change) Negative Binomial"
         ,y="log2(fold change) Gaussian/Poisson"
         ,title="Comparison of Fold Change Estimates") 
  # (negative binomial, ", variable,")
  pa <- 
  pa +
    geom_text(data=annot, aes(x=annotx
                              ,y=annoty
                              ,label=paste0("corr (",variable,"): "
                                     ,round(rhohat,2))
                              ,hjust=hjustvar
                              ,vjust=vjustvar
                              )
                              ,inherit.aes=FALSE) + 
    theme(
      legend.justification = 'left', 
      legend.position = 'bottom', legend.box = 'vertical', 
      legend.box.just = 'left') + 
    theme(legend.margin = margin(-5, 0, 0, 0))
  
  ####  Fold change comparisons
  palist[[ctt]] <- pa

  #### Signed -log10 pvalues  
  ovrplot[,sgnl10p:=-1*sign(log(fold_change)) * (-log10(p.value))]
  
  dm_p2 <- cast_and_melt(
    ovrplot
    ,castvar = c("modeltype")
    ,valuevar="sgnl10p"
    ,idvars = c("contrast", "cell_type_label", "frmla", "slide_use"
                ,"target", "imm_anno", "maxratio","cell_type_n", "ncells_1", "counts_1"
                , "cpc_1", "cpc_cell_type", "cell_type_cts")
    ,reflevel = "negative binomial RE"
  )
  pd2 <- dm_p2[!is.na(ncells_1)][cell_type_cts > 0][!is.na(value)]
  pd2[,cpc_cut:=cut(cpc_cell_type, c(0, 0.01, 0.05, 0.1, 100)
                   ,labels = c("< 0.01", "0.01-0.05", "0.05-0.1", "> 0.1")
                   , include.lowest=TRUE)]
  
  pd2[,lbl:=paste0(cell_type_label, " (", scales::comma(cell_type_n), " cells)")]
  pd2[,lbl:=factor(lbl, levels=pd2[order(cell_type_n)][,unique(lbl)])]
  fam_cols_d2 <- c("#1B9E77", "#D95F02", "#7570B3")
  names(fam_cols_d2) <- c("negative binomial RE","poisson RE", "gaussian RE")
  
  zthresh <- -qnorm(0.05/960/2)
  pthresh <- -log10(0.05/960)
  
  pp <- 
    ggplot(pd2[maxratio < 1][cpc_1 > 0.01][counts_1 > 10]#[cell_type_n > 5000]#[!is.infinite(value) ][! is.na(value)][!is.infinite(ranef_tissue_FALSE_nbinom2)]#[groupadj=="fixef_tissue"]
         ,aes(x=`negative binomial RE` #ranef_tissue_FALSE_nbinom2
              ,y =value 
              ,color=variable
              ,group=variable
              ,shape=cpc_cut
         )) + 
    geom_point(size=1) + 
    geom_abline(slope=1,intercept=0,lty=2) + 
    geom_smooth(method='lm',se=FALSE,alpha=0.5)  + 
    geom_hline(yintercept = pthresh,lty=2) + 
    geom_hline(yintercept = -pthresh,lty=2) + 
    geom_vline(xintercept = pthresh,lty=2) + 
    geom_vline(xintercept = -pthresh,lty=2) + 
    theme_bw() + 
    scale_color_manual(values=fam_cols_d2, guide=guide_legend(override.aes=list(size=4))
                       ,name="alternative model") + 
    scale_shape_manual(values=c(4,2,1, 7)
                       ,guide=guide_legend(override.aes=list(size=4))
                       ,name="counts per cell") + 
    coord_fixed() + 
    labs(x = " (effect signed) -log10(p) (negative binomial RE)"
         ,y = " (effect signed) -log10(p) (alternative model)"
         ,title="Comparison of effect sign and -log10(p)") 
  
  sigcl <- brewer.pal(8,"Dark2")[c(4,6)]
  
  txtd <-   
  pd2[maxratio < 1][cpc_1 > 0.01][counts_1 > 10][
      ((abs(`negative binomial RE`) > pthresh) & (abs(value) < pthresh)) | 
      ((abs(`negative binomial RE`) < pthresh) & (abs(value) > pthresh))
      ][abs(value) < 20][abs(`negative binomial RE`) < 20]
    
  txtd[,rnknb:=rank(-abs(`negative binomial RE`))]
  txtd[,rnkalt:=rank(-abs(value))]
  pp2 <- pp + 
         geom_label_repel(data=txtd[rnknb < 10 | rnkalt < 10][,head(.SD,1),by=target]
                         ,aes(label=paste0(target))
                         ,size=rplsize
                         ,hjust=0
                         ,show.legend = FALSE)
  
  
  pb <- 
  pp2 + 
    annotate("rect", xmin = pthresh, xmax = 25, ymin = -pthresh, ymax = pthresh,
             alpha = .1,fill = sigcl[1]) + 
    annotate("rect", xmin = -25, xmax = -pthresh, ymin = -pthresh, ymax = pthresh,
             alpha = .1,fill = sigcl[1]) + 
    annotate("text", x = pthresh+2, y = -pthresh + 0.5
             ,label='negative binomial\nsignificant, alternative not'
             ,size=3.5
             ,hjust=0,alpha = .8,color = sigcl[1])  +
    annotate("rect", xmin = -pthresh, xmax = pthresh, ymin = pthresh, ymax = 25,
             alpha = .1,fill = sigcl[2]) + 
    annotate("rect", xmin = -pthresh, xmax = pthresh, ymin = -25, ymax = -pthresh,
             alpha = .1,fill = sigcl[2]) + 
    annotate("text", x = -pthresh+1, y = -18,label='alternative model significant /\nnegative binomial not'
             ,size=3.5
             ,hjust=0,alpha = .8,color = sigcl[2]) 
  
    
  pblist[[ctt]] <- pb

  pdvenn <- pd2[maxratio < 1][cpc_1 > 0.01][counts_1 > 10]
  pdvenn[,.N,by=.(contrast, cell_type_label)]
  pdv <- copy(pdvenn)
  siglist <- list(
    "negative binomial RE" = pdv[abs(`negative binomial RE`) > zthresh][,.(contrast, target, cell_type_label)][,unique(.SD)][,paste0(target, ".", cell_type_label, ".", contrast)]
    ,"poisson RE" =       pdv[variable=="poisson RE"][abs(value) > zthresh][,paste0(target, ".", cell_type_label, ".", contrast)]
    ,"gaussian RE" =      pdv[variable=="gaussian RE"][abs(value) > zthresh][,paste0(target, ".", cell_type_label, ".", contrast)]
  )  
    ggv <- ggvenn(siglist
                  ,fill_color = unname(c(fam_cols_d2))
                  ,stroke_size=0.5
                  ,set_name_size = 4
                  ) + 
      labs(title = "Significant gene sets across contrasts\nby model type")
  

  pclist[[ctt]] <- ggv
  
}

pdall <- rbindlist(pdl)
annotall <- rbindlist(annotl)
annotall
annotall[,annotx2:=-Inf]
annotall[,annoty2:=-Inf]
pa_all <- 
ggplot(pdall[maxratio < 1][cpc_cell_type > 0.01][counts_1 > 10][value!=0]#[cpc_cell_type > 0.01]#[cts > 10]#[groupadj=="fixef_tissue"]
       ,aes(x=log2(nbinom2)
            ,y = log2(value)
            ,color=variable
            ,group=variable
            ,shape=cpc_cut
            # ,shape=contrast
       )) + 
  geom_point(size=1) + 
  theme_bw() + 
  scale_color_manual(values=fam_cols, guide=guide_legend(override.aes=list(size=4))
                     ,name="Alternative Distribution") + 
  scale_shape_manual(values=c(4,2,1, 7)
                     ,guide=guide_legend(override.aes=list(size=4))
                     ,name="counts per cell") + 
  geom_abline(slope=1,intercept=0) + 
  geom_smooth(method='lm',se=FALSE)  + 
  labs(x = "log2(fold change) Negative Binomial"
       ,y="log2(fold change) Gaussian/Poisson"
       ,title="Comparison of Fold Change Estimates")  + 
  facet_wrap(~cell_type_label)

annotall
annotall[,annotx2:=-Inf]
annotall[,annoty2:=-Inf]
annotall[,vjustvar2:=-1*vjustvar]
annotall[,hjustvar2:=hjustvar + 0.5]

pa_all2 <- 
pa_all +
  geom_text(data=annotall[,cpc_cut:=NA]
            , aes(
                   x=annotx2
                  ,y=annoty2
                  ,label=paste0("corr (",variable,"): "
                         ,round(rhohat,2))
                  ,hjust=hjustvar2
                  ,vjust=vjustvar2
                  )
            ) + 
  coord_fixed() + 
  theme(
    legend.justification = 'left', 
    legend.position = 'bottom', legend.box = 'vertical', 
    legend.box.just = 'left') + 
  theme(legend.margin = margin(-5, 0, 0, 0)) 

pa_all2
#################################################################################
#################################################################################
###             Supplementary Figure 1               ###########################
#################################################################################
#################################################################################

print(pa_all2)
ggsave("plots/fixedsize__supl_effectsize_by_family_allcelltypes.png"
       ,width=823*sqrt(30)
       ,height=1120*sqrt(30)
       ,dpi=450
       ,units='px')



################################################################################
### Suppl Figure 2, all celltypes signed -log10(p) b/w families
################################################################################
pb_legend <- cowplot::get_legend(pblist$`T CD4` + coord_fixed() + 
                                #   theme(legend.box.margin = margin(0, 0, 0, 0))
                                   theme(legend.position="right")
                                 )
pblist_comb <- 
lapply(names(palist), function(ct){
  pblist[[ct]] + 
    coord_fixed(xlim=c(-20,20), ylim=c(-20, 20)) + 
    theme(legend.position="none") + 
    labs(title="",subtitle=paste0(ct)) + 
    theme(plot.margin = unit(c(0.5, 0.5, 0.5, 0.5), "cm"))
})
names(pblist_comb) <- names(pblist)
pblist_comb[["legend"]] <- pb_legend

cp <- cowplot::plot_grid(plotlist= pblist_comb, nrow = 3)  + 
  theme(plot.background = element_rect(fill = 'white',color='white')) + 
  cowplot::panel_border(remove=TRUE)
print(cp)


ggsave("plots/supl_signedl10p_by_family_allcelltypes.png"
       ,width=1109*sqrt(20)
       ,height=1049*sqrt(20)
       ,dpi=350
       ,units='px')

ncells_lbls <- c('Treg: (4,361 cells)'
                 ,'T CD8: (9,433 cells)'
                 ,'plasmablast: (21,557 cells)'
                 ,'T CD4: (29,905 cells)'
                 ,'macrophage: (34,590 cells)'
                 ,'neutrophil: (42,219 cells)'
                 ,'tumor: (227,812 cells)'
                 )
names(ncells_lbls) <- gsub(":.*$", "", ncells_lbls)

pclist$`T CD4`
pclist_comb <- 
lapply(names(ncells_lbls), function(ct){
  pclist[[ct]] + 
  labs(title="",subtitle=paste0(ncells_lbls[ct])) + 
  theme(plot.margin = unit(c(0.5, -100, 0.5, -100), "cm"))
})
pclist_comb[[1]]
cpv <- cowplot::plot_grid(plotlist= pclist_comb, nrow = 3)  + 
  theme(plot.background = element_rect(fill = 'white',color='white')) + 
  cowplot::panel_border(remove=TRUE)
print(cpv)
ggsave("plots/supl_vennsig_by_family_allcelltypes.png"
       ,width=1097*sqrt(20)
       ,height=1167*sqrt(20)
       ,dpi=350
       ,units='px')

ncells_lbls

##########################################
#### Computation time plot ##############
##########################################
ovrplotuse <- fread("nsclc_de_results/ovr.plotting.csv")
pd <- copy(ovrplotuse) #fread("plotdata/ovr.comparefamilies.clean.csv")
pd <- pd[,ncells_lbl:=paste0(cell_type_label, ": (", scales::comma(cell_type_n), " cells)")]
pd[,ncells_lbl:=factor(ncells_lbl, levels=pd[order(cell_type_n),unique(ncells_lbl)])]
pd[,variablef:=paste0(groupadj, "_FALSE_", fam)]
pd[variablef=="ranef_tissue_FALSE_poisson",   variablef:="poisson RE"]
pd[variablef=="ranef_tissue_FALSE_gaussian",  variablef:="gaussian RE"]
pd[variablef=="ranef_tissue_FALSE_nbinom2",   variablef:="negative binomial RE"]
pd[,groupadj_lbl:=paste0("random effect per sample")]
fam_cols_d3 <- c(fam_cols_d2, brewer.pal(8,"Dark2")[8])
names(fam_cols_d3)[length(fam_cols_d3)] <- "negative binomial RE"

comptime_plot <- 
  ggplot(pd
         ,aes(x=ncells_lbl,y=elapsed,fill = variablef)) + 
  geom_boxplot(outlier.size=0.2) + coord_flip() + 
  theme_bw() + 
  theme(axis.text.x=element_text(face='bold')) + 
  theme(axis.text.y=element_text(face='bold')) + 
  scale_y_continuous(trans='log10',name = "elapsed (seconds)",labels=scales::comma
  )  + 
  scale_fill_manual(values=fam_cols_d2, guide=guide_legend(override.aes=list(size=4))
                     ,name="model type") + 
  labs(title='Model fitting time per gene'
       ,x="cell type (# of cells)")  + 
  theme(legend.position="bottom") + 
  theme(axis.text.y=element_text(angle=45))

comptime_plot



##########################################
###  Figure 2      #######################
##########################################
addcomptime <- FALSE
addcomptime <- TRUE
overallplot <- list()
ctt <- "macrophage"
pab <-   
cowplot::plot_grid(palist[[ctt]]  + 
                     theme(legend.position="none") + 
                     theme(text=element_text(face='bold',size=12)) + 
                     coord_fixed()
                  # ,NULL
                    , pblist[[ctt]] + coord_fixed(xlim=c(-20,20), ylim=c(-20,20)) + 
                     theme(text=element_text(face='bold',size=12))
                    ,rel_widths = c(1, 1.05)
                   ,labels = 'auto', label_size=24) 

if(addcomptime){
  pcd <-   
  cowplot::plot_grid(pclist[[ctt]] + 
                       theme_bw() +
                       theme(axis.line = element_blank(),
                             axis.text.x=element_blank(),
                             axis.text.y=element_blank(),
                             axis.ticks.x=element_blank(),
                             axis.ticks.y=element_blank(),
                             axis.title.x=element_blank(),
                             axis.title.y=element_blank(),
                             panel.grid.major = element_blank(),
                             panel.grid.minor = element_blank(),
                             panel.border = element_blank(),
                             panel.background = element_blank()
                             ,text=element_text(size=12, face='bold')
                             ) + 
                       theme(plot.margin = unit(c(0,0,0,0), "cm"))
                     ,comptime_plot + theme(plot.margin = unit(c(0,0.2,0,0), "cm")) + 
                      theme(text=element_text(size=12, face='bold'))
                     ,rel_widths = c(0.7,1)
                     ,labels = c('c','d'), label_size = 24) 
  
} else {
  pcd <-   
  cowplot::plot_grid(pclist[[ctt]] + 
                       theme_bw() +
                       theme(axis.line = element_blank(),
                             axis.text.x=element_blank(),
                             axis.text.y=element_blank(),
                             axis.ticks.x=element_blank(),
                             axis.ticks.y=element_blank(),
                             axis.title.x=element_blank(),
                             axis.title.y=element_blank(),
                             panel.grid.major = element_blank(),
                             panel.grid.minor = element_blank(),
                             panel.border = element_blank(),
                             panel.background = element_blank()) + 
                       theme(plot.margin = unit(c(0,0,0,0), "cm"))
                     ,NULL
                     ,rel_widths = c(0.7,1)
                     ,labels = c('c','')) 
  
}

pabcd <-   
cowplot::plot_grid(pab 
                   ,pcd + theme(plot.margin = unit(c(0,0,0,0), "cm"))
                   ,nrow = 2
                   ) 
title_gg <- ggplot() + theme_minimal() + 
  labs(title = paste0("Comparison between parametric distributions in DE regression models; ", paste0("(", ctt, " cell type)"))
       )  + theme(plot.title=element_text(hjust=0.5
                                     )
                  ,plot.subtitle = element_text(hjust=0.4)
                  )
overallplot[[ctt]] <- 
cowplot::plot_grid(pabcd,nrow=1)

print(
  overallplot[["macrophage"]] + 
  theme(plot.background = element_rect(fill = 'white',color='white')) + 
  cowplot::panel_border(remove=TRUE)
  )

ggsave("plots/Figure2.png"
       ,width=1179*sqrt(20)
       ,height=776*sqrt(20)
       ,dpi=350
       ,units='px')

