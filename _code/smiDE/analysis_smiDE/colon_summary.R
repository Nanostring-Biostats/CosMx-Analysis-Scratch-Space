# devtools::install_github("Winnie09/Palo")
library(data.table)
library(ggplot2)
library(RColorBrewer)
library(ggrepel)

orm_metrics <- fread("colon_de_results/orm_metrics_colon.csv")
load("colon_data/coloncancer.RData")

RhpcBLASctl::blas_set_num_threads(1)
metainfo <- data.table(cbind(annot, clust))

#targset <- as.numeric(args[1])

### pre de
pre_de_bcell <- 
  smiDE::pre_de(counts = Matrix::t(raw)
                ,normalized_data = Matrix::t(raw)
                ,weight_colname = NULL
                ,metadata = metainfo
                ,ref_celltype = "B-cell"
                ,cell_type_metadata_colname = "clust"
                ,sdimx_colname = "x_slide_mm"
                ,sdimy_colname = "y_slide_mm"
                ,mm_radius = 0.015 ### 15 microns
                ,aggregation = "sum"
                ,adjacencies_only = FALSE)

#pre_de_bcell$nblist$adjacency_counts_by_ct[1:10] ### count # of cells of each 
#                                                 ### neighbor cell type

bcell_cellid <- metainfo[clust=="B-cell", cell_ID]
pre_de_bcell$nblist$adjacency_counts_by_ct[match(bcell_cellid
                                                 ,cell_ID)
                                           ,summary(rowSums(.SD))
                                           ,.SDcols=setdiff(colnames(pre_de_bcell$nblist$adjacency_counts_by_ct), "cell_ID")]

immune_types <- c("B-cell", "macrophage", "mDC", "neutrophil"
                  ,"pDC", "mast", "monocyte", "NK", "plasmablast"
                  ,"T CD8 memory", "plasma.cell", "plasmablast"
                  ,"T CD4 naive", "Treg")

setdiff(metainfo[,unique(clust)], immune_types)

celltype_nbrs <- 
  copy(pre_de_bcell$nblist$adjacency_counts_by_ct)[,immune_ct_neighbors:=rowSums(.SD)
                                                   ,.SDcols=c(immune_types)
  ]

metainfo_bcell <- merge(metainfo[clust=="B-cell"]
                        ,celltype_nbrs[,.(cell_ID, immune_ct_neighbors)]
                        , by="cell_ID")

ggplot(metainfo, aes(x_slide_mm, y_slide_mm, color=factor(fov))) + 
  geom_point(size=0.1) + 
  geom_point(data=metainfo[clust=="B-cell"], size=0.1) + 
  coord_fixed() + 
  theme_bw() + 
  scale_color_manual(values=rep(unname(pals::alphabet()),3)
#  scale_color_manual(values=ctpal[levs]
 #                    ,breaks=levs
                     ,guide=guide_legend(override.aes=list(size=4))) + 
  geom_text(data=metainfo[,lapply(.SD, mean),.SDcols=c("x_slide_mm", "y_slide_mm"), by=fov]
            ,aes(x=x_slide_mm, y=y_slide_mm, label=fov), inherit.aes = FALSE, color='black') 


metainfo[,tissgrp:="ll"]
metainfo[fov >=60, tissgrp:="ur"]
metainfo[,bcell_nbr:=0]
metainfo[match(pre_de_bcell$cell_adjacency_dt[from %in% metainfo_bcell[["cell_ID"]] | 
                                                to %in% metainfo_bcell[["cell_ID"]]
                                              ,unique(c(from, to))]
               ,cell_ID), bcell_nbr:=1]
metainfo[,.N,by=.(bcell_nbr)]

palcols <- unique(
            c(brewer.pal(12, "Paired")
              ,brewer.pal(8, "Dark2")
              ,brewer.pal(8, "Accent")
              ,brewer.pal(9,"Set3")[2:6])
          )

palcols <- palcols[1:length(unique(metainfo[["clust"]]))]
allcols_discrete <- c(
  "#7FC97F"  ,"#BEAED4" ,"#FDC086" ,"#FFFF99"
  ,"#386CB0" ,"#F0027F" ,"#BF5B17" ,"#666666"
  ,"#1B9E77" ,"#D95F02" ,"#7570B3" ,"#E7298A"
  ,"#66A61E" ,"#E6AB02" ,"#A6761D" ,"#666666"
  ,"#A6CEE3" ,"#1F78B4" ,"#B2DF8A" ,"#33A02C"
  ,"#FB9A99" ,"#E31A1C" ,"#FDBF6F" ,"#FF7F00"
  ,"#CAB2D6" ,"#6A3D9A" ,"#FFFF99" ,"#B15928"
  ,"#FBB4AE" ,"#B3CDE3" ,"#CCEBC5" ,"#DECBE4"
  ,"#FED9A6" ,"#FFFFCC" ,"#E5D8BD" ,"#FDDAEC"
  ,"#F2F2F2" ,"#B3E2CD" ,"#FDCDAC" ,"#CBD5E8"
  ,"#F4CAE4" ,"#E6F5C9" ,"#FFF2AE" ,"#F1E2CC"
  ,"#CCCCCC" ,"#E41A1C" ,"#377EB8" ,"#4DAF4A"
  ,"#984EA3" ,"#FF7F00" ,"#FFFF33" ,"#A65628"
  ,"#F781BF" ,"#999999" ,"#66C2A5" ,"#FC8D62"
  ,"#8DA0CB" ,"#E78AC3" ,"#A6D854" ,"#FFD92F"
  ,"#E5C494" ,"#B3B3B3" ,"#8DD3C7" ,"#FFFFB3"
  ,"#BEBADA" ,"#FB8072" ,"#80B1D3" ,"#FDB462"
  ,"#B3DE69" ,"#FCCDE5" ,"#D9D9D9" ,"#BC80BD"
  ,"#CCEBC5" ,"#FFED6F"
)
set.seed(240404)
ctpal <- 
Palo::Palo(metainfo[,.(x_slide_mm, y_slide_mm)]
           ,cluster=metainfo$clust
           ,palette = unique(allcols_discrete)[1:length(unique(metainfo$clust))]
           )
ctpal["B-cell"] <- "#FFFF99"
ctpal["epithelial.normal.crypts"] <- "#E6AB02"
levs <- metainfo[,.N,by=.(clust)][order(-N),clust]
levs <- c("B-cell", setdiff(levs, "B-cell"))
metainfo[,cell_type:=factor(clust, levels=levs)]
plist <- 
lapply(split(metainfo, by="tissgrp"), function(xx){
  ggplot(xx[clust!="B-cell"], aes(x_slide_mm, y_slide_mm, color=cell_type)) + 
    geom_point(size=0.1) + 
    geom_point(data=xx[clust=="B-cell"], size=0.1) + 
    coord_fixed() + 
    theme_bw() + 
    #scale_color_manual(values=rep(unname(pals::alphabet()),2)
    scale_color_manual(values=ctpal[levs]
                       ,breaks=levs
                       ,guide=guide_legend(override.aes=list(size=4)))
    
})
cowplot::plot_grid(plist[[1]], plist[[2]] + theme(legend.position = "none")
                   ,rel_widths = c(8,4))

metainfo[,ct:="other"]
metainfo[clust %in% immune_types,ct:="other immune cell type"]
metainfo[clust == "B-cell",ct:="B cell"]
metainfo[,ct:=factor(ct, levels=c("other", "other immune cell type", "B cell"))]
ctpal2 <- c("grey50", "#33A02C", "black")
names(ctpal2) <- c("other", "other immune cell type", "B cell")


plist <- 
lapply(split(metainfo, by="tissgrp"), function(xx){
  ggplot(xx[clust!="B-cell"][order(ct)], aes(x_slide_mm, y_slide_mm, color=ct)) + 
    geom_point(size=0.05, alpha=0.5) + 
    geom_point(data=xx[clust=="B-cell"], size=0.3) + 
    coord_fixed() + 
    theme_bw() + 
    scale_color_manual(values=rev(ctpal2)
                       ,breaks=rev(names(ctpal2))
                       ,name="cell type"
                       ,guide=guide_legend(override.aes=list(size=4, alpha=1))) + 
    theme(legend.position="bottom")
    
})
cowplot::plot_grid(plist[[1]] + theme(plot.margin = unit(c(0.1,-2,0,0), "cm")) + theme(legend.position="none")
                   ,plist[[2]] + theme(plot.margin = unit(c(0,0.1,0,-2), "cm")) + 
                     theme(legend.position = "none")
                   ,rel_widths = c(8,4))


pll <- plist[[1]] + theme(plot.margin = unit(c(0.1,-2,0,0), "cm")) + theme(legend.position="none")
pll2 <- plist[[2]] + theme(plot.margin = unit(c(0,0.1,0,-2), "cm")) + theme(legend.position = "none")
pcelltype <- 
cowplot::plot_grid(pll + theme(plot.margin = unit(c(0,-2,0,0), "cm"))
                   ,cowplot::plot_grid(pll2,NULL, rel_heights = c(3, 4.2), nrow=2) + 
                      theme(plot.margin = unit(c(0,0.1,0,-2), "cm"))
                   ,rel_widths = c(3, 2.5))

pcelltype

metainfo_bcell[,.N,by=.(immune_ct_neighbors)]
metainfo_bcell[,imm_cut:=cut(immune_ct_neighbors
                       ,quantile(immune_ct_neighbors, seq(0,1,1/7))
                       ,include.lowest=TRUE)]

metainfo_bcell[,.N,by=.(imm_cut)]
metainfo_bcell[,tissgrp:="ll"]
metainfo_bcell[fov >=60, tissgrp:="ur"]
cutpal <- colorRampPalette(brewer.pal(11, "RdYlBu"))(7)
names(cutpal) <- metainfo_bcell[,levels(imm_cut)]
plist_immcut <- 
lapply(split(metainfo_bcell, by="tissgrp"), function(xx){
  ggplot(xx[order(imm_cut)]
         ,aes(x_slide_mm, y_slide_mm, fill=imm_cut)) + 
    geom_point(size=2
               ,pch=21
               ,alpha=1
               ) + 
    coord_fixed() + 
    theme_bw() + 
    scale_fill_manual(name = "# of immune cell neighbors"
                      ,values=cutpal, guide=guide_legend(
                                                         override.aes=list(size=4,alpha=1))) + 
    theme(legend.position="right")
})

rightside3 <- cowplot::plot_grid(plist_immcut$ll + theme(legend.position="none") +
                                   theme(axis.text = element_blank(),axis.ticks = element_blank(), axis.title = element_blank()) + 
                                         theme(plot.margin = unit(c(0.1,0,0,0), "cm"))
                                       ,cowplot::plot_grid(plist_immcut$ur + theme(legend.position="none") +
                                         theme(axis.text = element_blank(), axis.ticks=element_blank() ,axis.title = element_blank()) + 
                                         theme(plot.margin = unit(c(0,0,0,0), "cm")), NULL, rel_widths = c(2.5,0.1))
                                       , nrow=2
                                       ,rel_heights = c(6, 3.3)
                                ,align = 'v'
          )
rightside3


leftside <- 
cowplot::plot_grid(cowplot::plot_grid(plist$ll + theme(legend.position="none") 
                                      ,plist$ur + theme(legend.position="none")
                                      ,nrow=2, rel_heights = c(6,3))
                   )

cowplot::plot_grid(leftside, rightside3, nrow=1)

leg1 <- cowplot::get_plot_component(plist[[1]] + theme(legend.position = "right"), 'guide-box-right', return_all=TRUE)
leg2 <- cowplot::get_plot_component(plist_immcut[[1]] , 'guide-box-right', return_all=TRUE)


combleg <- cowplot::plot_grid(cowplot::ggdraw(leg1) + theme(plot.margin = unit(c(0,-5,0,-5), "cm"))
                              ,cowplot::ggdraw(leg2)+ theme(plot.margin = unit(c(0,-5,0,-5), "cm"))
                              , nrow=2)

plot_ab <- 
cowplot::plot_grid(leftside
                   ,combleg
                   ,rightside3
                   ,nrow=1
                   ,labels = c("a", "", "b")
                   ,rel_widths = c(1, 0.4, 1)) + 
  theme(plot.background = element_rect(fill = 'white',color='white')) + 
  cowplot::panel_border(remove=TRUE)
print(plot_ab)
ggsave("colon_data/deplot1.png"
       ,width=938*sqrt(25)
       ,height=904*sqrt(25)
       ,dpi=450
       ,units = "px")

dev2f <- list.files("colon_de_results", pattern="de_seminaive",full.names=TRUE)
dev1f <- list.files("colon_de_results", pattern="de_naive",full.names=TRUE)
dev3f <- list.files("colon_de_results", pattern="de_spatial",full.names=TRUE)

deobj <- readRDS(dev2f[1])
for(ff in 2:length(dev2f)){
  message(ff / length(ff))
  deobjff <- readRDS(dev2f[ff])
  deobj$results <- c(deobj$results, deobjff$results)
  deobj$targets <- c(deobj$target, deobjff$targets)
}
saveRDS(deobj, file="colon_data/seminaive_de_results.rds")

deobj <- readRDS(dev1f[1])
for(ff in 2:length(dev1f)){
  message(ff / length(ff))
  deobjff <- readRDS(dev1f[ff])
  deobj$results <- c(deobj$results, deobjff$results)
  deobj$targets <- c(deobj$target, deobjff$targets)
}
saveRDS(deobj, file="colon_data/naive_de_results.rds")

deobj <- readRDS(dev3f[1])
for(ff in 2:length(dev3f)){
  message(ff / length(ff))
  deobjff <- readRDS(dev3f[ff])
  deobj$results <- c(deobj$results, deobjff$results)
  deobj$targets <- c(deobj$target, deobjff$targets)
}
saveRDS(deobj, file="colon_data/inlaspatialaware_de_results.rds")


ms2 <- rbindlist(msl)
pw2 <- rbindlist(pwl)
emm2 <- rbindlist(emml)


dev2f <- list.files("colon_data/de_results", pattern="de_v2",full.names=TRUE)
dev1f <- list.files("colon_data/de_results", pattern="naive",full.names=TRUE)
dev3f <- list.files("colon_data/de_results", pattern="de_v3",full.names=TRUE)

pwl <- emml <- msl <- list()
for(ff in 1:length(dev2f)){
  deobj <- readRDS(dev2f[ff])
  pwl[[ff]] <- smiDE::results(deobj, "pairwise")[[1]]
  emml[[ff]] <- smiDE::results(deobj, "emmeans")[[1]]
  msl[[ff]] <- smiDE::results(deobj, "model_summary")[[1]]
}
ms2 <- rbindlist(msl)
pw2 <- rbindlist(pwl)
emm2 <- rbindlist(emml)

pwl <- emml <- msl <- list()
for(ff in 1:length(dev3f)){
  deobj <- readRDS(dev3f[ff])
  pwl[[ff]] <- smiDE::results(deobj, "pairwise")[[1]]
  emml[[ff]] <- smiDE::results(deobj, "emmeans")[[1]]
  msl[[ff]] <- smiDE::results(deobj, "model_summary")[[1]]
}
ms3 <- rbindlist(msl)
pw3 <- rbindlist(pwl)
emm3 <- rbindlist(emml)


pwl <- emml <- msl <- list()
for(ff in 1:length(dev1f)){
  deobj <- readRDS(dev1f[ff])
  pwl[[ff]] <- smiDE::results(deobj, "pairwise")[[1]]
  emml[[ff]] <- smiDE::results(deobj, "emmeans")[[1]]
  msl[[ff]] <- smiDE::results(deobj, "model_summary")[[1]]
}


ms1 <- rbindlist(msl)
pw1 <- rbindlist(pwl)
emm1 <- rbindlist(emml)


srel <- list()
for(ff in 1:length(dev3f)){
  message(ff)
  deobj <- readRDS(dev3f[ff])
  srel[[ff]] <- smiDE::results(deobj, "spatial_random_effect")[[1]]
}
sre <- rbindlist(srel)


vnaive <- 
smiDE::volcano(pw1[term=="immune_ct_neighbors"][!grepl("Custom", target)]
               ,fold_change_cutoff = 2**0.5
               , interactive=FALSE
               ,limit_labels = TRUE)[[1]]

vnaive <- 
vnaive + labs(title="Naive DE Analysis"
              ,subtitle="NO Contamination gene filtering,\nNO covariate adjustment using gene expression in neighboring cells,\nNO spatial random effect"
              ,x=paste0("# of immune neighbors (9.458)", 
    " ", sprintf("←"), " log", sprintf("₂"), "(FC) ", 
    sprintf("→"), " ", "# of immune neighbors (14.54)"))
vnaive

v2 <- 
smiDE::volcano(pw2[term=="immune_ct_neighbors"][!grepl("Custom", target)]
               ,fold_change_cutoff = 2**0.5
               , interactive=FALSE
               ,limit_labels = TRUE)[[1]]

v2 <- 
v2 + labs(title="Semi-Naive DE Analysis"
              ,subtitle="With Contamination gene filtering,\nWith Covariate adjustment using gene expression in neighboring cells,\nNO spatial random effect"
              ,x=paste0("# of immune neighbors (9.458)", 
    " ", sprintf("←"), " log", sprintf("₂"), "(FC) ", 
    sprintf("→"), " ", "# of immune neighbors (14.54)"))

v2 + theme(plot.title = element_text(size=8), plot.subtitle=element_text(size=8))


pw1p <- 
merge(pw1[term=="immune_ct_neighbors"][!grepl("Custom", target)]
      ,orm_metrics[clust=="B-cell",.(target, contam_ratio=ratio)]
      ,by="target"
      )

pw2p <- 
merge(pw1p, pw2[term=="immune_ct_neighbors"], by=c("target"), suffixes=c("_v1", "_v2"))

pw3p <- merge(pw2p, pw3[term=="immune_ct_neighbors"], by="target")

pw3p <- 
merge(pw3p, ms3[term=="Stdev for s"][,.(spatial_re_sd_mean = mean,target)]
      ,by="target")

pw3p[,lcl_bonf:=`9.255831e-06quant`]
pw3p[,ucl_bonf:=`0.9999907quant`]
pw3p[,lcl:=`0.025quant`]
pw3p[,ucl:=`0.975quant`]

pw3p[,sig_bonf:=as.numeric(sign(log(lcl_bonf))==sign(log(ucl_bonf)))]
pw3p[,sig:=as.numeric(sign(log(lcl))==sign(log(ucl)))]
pw3p[,showbound:=`0.5quant`]
pw3p[sig==1 & `0.975quant` > 1,showbound:=`0.025quant`]
pw3p[sig==1 & `0.975quant` < 1,showbound:=`0.975quant`]
pw3p[`0.975quant` < 1, showbound:=`0.975quant`]
pw3p[,sigcat:="nonsig"]
pw3p[sig==1,sigcat:="significant at 95% CI"]
pw3p[sig_bonf==1,sigcat:="significant at 95% Bonferroni adjusted CI"]

pw3p[,.(sig*abs(`0.5quant`-1), `0.5quant`, log2(`0.5quant`), sig*abs(log2(`0.5quant`)-0))]
pw3p[,yval1:=sig*abs(log2(`0.5quant`)-0)]
pw3p[,yval2:=sig*abs(`0.5quant`-1)]
pw3p[,yval3:=sig*abs(`0.5quant`*(14.54-9.458)-1)]

sigpal <- c("grey", brewer.pal(9,"Pastel1")[2],scales::muted("red"))
names(sigpal) <- c("nonsig", "significant at 95% CI", "significant at 95% Bonferroni adjusted CI")
pw3p[,version:="spatial-aware (INLA)"]
vinla <- 
ggplot(pw3p, aes(x=log2(showbound)
                 ,y = yval2
                 ,color=sigcat)) + 
  geom_point() + 
  theme_bw(base_size=20) + 
  scale_color_manual(values=sigpal
                       ,guide=guide_legend(override.aes=list(size=4),nrow=2)
                     ,breaks=names(sigpal)[2:3]
                     ,name="Significance") + 
  geom_label_repel(data=pw3p[sigcat=="significant at 95% Bonferroni adjusted CI"]
                   ,aes(label=target)
                   ,show.legend=FALSE
                   ) + 
  labs(x="log2 (Fold Change Credible Interval Bound closest to 1)"
       ,y="Fold Change (+ Positive, - Negative)") + 
  scale_y_continuous(breaks=seq(0,0.8, 0.1)
                     ,labels=paste0("+ ", c(1 + seq(0,0.8, 0.1)), "\n-", c(1 - seq(0,0.8, 0.1)))) + 
  facet_wrap(~version, labeller=label_both)

vinla

pw12 <- 
rbindlist(list(pw1[,version:="naive"]
               ,pw2[,version:="semi-naive"])
          ,use.names=TRUE,fill=TRUE)

v12 <- 
smiDE::volcano(pw12[term=="immune_ct_neighbors"][!grepl("Custom", target)]
               ,fold_change_cutoff = 2**0.5
               , interactive=FALSE
               ,limit_labels = TRUE)[[1]] + 
  facet_wrap(~version,labeller=label_both)

v12 + labs()
g1 = ggplotGrob(v12)
gtable_show_layout(g1)


xaxis <- g1[c(12),]
ggplotify::as.ggplot(xaxis)

g1_am0 <- g1[,c(5:9)]
g1p <- ggplotify::as.ggplot(g1_am0)
g2_am0 <- g1[,c(9:12)]

g2p <- ggplotify::as.ggplot(g2_am0)

gab <- 
cowplot::plot_grid(g1p  #+ labs(title="Naive model with no control variable adjustment or filtering") + theme(plot.title=element_text(hjust=0.5, vjust=-8, size=20))
                   ,g2p #+ labs(title="Semi-Naive")
                   ,rel_widths = c(1,0.8)
                   ,labels=c("a", "b")
                   ,label_size = 20)
gab
gabp <- 
cowplot::plot_grid(gab + theme(plot.margin = unit(c(0,0,-1,0), "cm"))
                   ,xaxis #+ theme(plot.margin = unit(c(0,0,1,0), "cm"))
                   , rel_heights = c(1, 0.05), nrow=2)



cowplot::plot_grid(gabp
                   ,cowplot::plot_grid(NULL
                                       ,vinla + theme(legend.position = "bottom")
                                       ,NULL
                                       ,rel_heights = c(0.2, 1, 0.1)
                                       ,nrow=3
                                       )
                   ,nrow=1, rel_widths = c(2, 0.8))

gabcp <- 
cowplot::plot_grid(gabp + 
  theme(plot.background = element_rect(fill = 'white',color='white')) + 
  cowplot::panel_border(remove=TRUE)
                   ,vinla + theme(legend.position = "top"
                                  ,axis.title.x = element_text(size=16)
                                  )
                   ,labels=c("", "c"), label_size = 20
                   ,nrow=1, rel_widths = c(2, 1))

gabcp + 
  theme(plot.background = element_rect(fill = 'white',color='white')) + 
  cowplot::panel_border(remove=TRUE)

ggsave("colon_data/deplot2.png"
       ,dpi=450
       ,width=7.54 * sqrt(80)#sqrt(100)
       ,height=700 * sqrt(80)#sqrt(100)
       ,units = "px"
       )
ggsave("colon_data/deplot22.png"
       ,dpi=450
       ,width=19*sqrt(15)
       ,height=8*2.54/3 *sqrt(15)
       ,units = "cm"
       )



pmiddle <- 
cowplot::plot_grid(g1p  #+ labs(title="Naive")
                   ,g2p #+ labs(title="Semi-Naive")
                   ,rel_widths = c(1,0.8)
                   ,labels=c("b", "c")
                   ,label_size = 20)

pmiddle <- 
cowplot::plot_grid(pmiddle + theme(plot.margin = unit(c(0,0,-1,0), "cm"))
                   ,xaxis #+ theme(plot.margin = unit(c(0,0,1,0), "cm"))
                   , rel_heights = c(1, 0.05), nrow=2)


#pmiddle <- 
#cowplot::plot_grid(pmiddle
#                   ,cowplot::plot_grid(NULL
#                                       ,vinla + theme(legend.position = "bottom")
#                                       ,NULL
#                                       ,rel_heights = c(0.2, 1, 0.1)
#                                       ,nrow=3
#                                       )
#                   ,nrow=1, rel_widths = c(2, 0.8))
#
pmiddle <- 
cowplot::plot_grid(pmiddle + 
  theme(plot.background = element_rect(fill = 'white',color='white')) + 
  cowplot::panel_border(remove=TRUE)
                   ,vinla + theme(legend.position = "top"
                                  ,axis.title.x = element_text(size=14))
                   ,labels=c("", "d"), label_size = 20
                   ,nrow=1, rel_widths = c(2, 1))

pmiddle <- pmiddle + 
  theme(plot.background = element_rect(fill = 'white',color='white')) + 
  cowplot::panel_border(remove=TRUE)





#
#
#
#
#
#cowplot::plot_grid(vnaive + scale_y_continuous(limits=ylim2), v2, vinla, nrow=1)
#
#ggplot2::layer_scales(vnaive)$y$range$range
#ylim2 <- ggplot2::layer_scales(v2)$y$range$range
#
#
#
#
#ggplot(pw2p, aes(x=-log10(p.value_v2) / -log10(p.value_v1),
#                 y=contam_ratio)) + 
#  geom_point()
#library(scales)
#ggplot(pw2p, aes(x=-log10(p.value_v2)
#                 ,y=-log10(p.value_v1),
#                 fill=contam_ratio)) + 
#  geom_point(pch=21,size=3) + 
#  theme_bw() + 
##  scale_color_gradient2() + 
#  scale_fill_gradient2(high=muted("red"), mid="white", low=muted("blue"), midpoint = median(pw2p[["contam_ratio"]])) + 
#  geom_abline(slope=1,intercept=0,lty=2) + 
#  scale_x_continuous(trans='log10', n.breaks=10) + 
#  scale_y_continuous(trans='log10', n.breaks=10)  
#
#ggplot(pw2p, aes(x=fold_change_v1
#                 ,y=fold_change_v2,
#                 fill=contam_ratio)) + 
#  geom_point(pch=21,size=3) + 
#  theme_bw() + 
##  scale_color_gradient2() + 
#  scale_fill_gradient2(high=muted("red"), mid="white", low=muted("blue"), midpoint = median(pw2p[["contam_ratio"]])) + 
#  geom_abline(slope=1,intercept=0,lty=2) + 
#  scale_x_continuous(trans='log2', n.breaks=10) + 
#  scale_y_continuous(trans='log2', n.breaks=10)  
  
#pw3p[,sig:=as.numeric(sign(log(`9.255831e-06quant`))==sign(log(`0.9999907quant`)))]

lbld <- rbindlist(list(pw3p[sig==1][order(`0.5quant`)][,head(.SD,10)]
                       ,pw3p[sig==1][order(`0.5quant`)][,tail(.SD,10)]
                       ,pw3p[sig_bonf==1]
                       ,pw3p[p.value_v2 < 0.05/.N][order(-abs(fold_change_v2))]
                       )
                  )[,unique(.SD)]

lbld <- rbindlist(list(pw3p[sig==1][order(`0.5quant`)][,head(.SD,10)]
                       ,pw3p[sig==1][order(`0.5quant`)][,tail(.SD,10)]
                       ,pw3p[sig_bonf==1]
                       ,pw3p[p.value_v2 < 0.05/.N][order(-abs(fold_change_v2))][1:10]
                       )
                  )[,unique(.SD)]
lbld[fold_change_v2 < `0.5quant`,ndg:=0.1 - fold_change_v2]
lbld[fold_change_v2 > `0.5quant`,ndg:=1.8 - fold_change_v2]


spatialp <- 
ggplot(pw3p, aes(x=fold_change_v2
                 ,y=`0.5quant`,
                 fill=spatial_re_sd_mean, shape=)) + 
  geom_point(data=pw3p[sig_bonf==1]
             ,aes(shape=0),size=6,stroke=1.2) + 
  geom_point(data=pw3p[p.value_v2 < 0.05/.N]
             ,aes(shape=8),size=3,stroke=1.2) + 
  geom_point(aes(shape=21),size=3) + 
  scale_shape_identity(guide=guide_legend(override.aes=list(size=6))
                       ,breaks=c(0,8)
                       ,name=""
                       ,labels=c("9 Significant genes\n(INLA spatial RE model)"
                                 ,"114 Significant genes\n(non-spatial semi-naive model)")) + 
  theme_bw(base_size = 20) + 
  labs(x="Fold change (without spatial RE)"
       ,y ="Fold change (INLA, with spatial RE)") + 
#  scale_color_gradient2() + 
  scale_fill_gradient2(name = "Std Dev of\nSpatial\nRandom Effect"
                       ,high=scales::muted("red"), mid="white", low=scales::muted("blue")
                       #,midpoint = median(pw3p[["spatial_re_sd_mean"]])
                       ,midpoint = mean(pw3p[["spatial_re_sd_mean"]])
                       ) + 
  geom_abline(slope=1,intercept=0,lty=2) + 
  scale_x_continuous( n.breaks=20) + 
  scale_y_continuous( n.breaks=20) + 
  #geom_label_repel(data=lbld[sig_bonf!=1]
  #                 ,aes(label=target)
  #                 ,max.overlaps=20
  #                 ,nudge_x = lbld[sig_bonf!=1,ndg]
  #                 ) + 
  geom_label_repel(data=lbld
                   ,aes(label=target)
                   ,max.overlaps=20
                   ,nudge_x = lbld[,ndg]
                   ,direction="y"
                   ) + 
#  geom_label_repel(data=lbld[sig_bonf==1]
#                   ,aes(label=target)
#                   ,max.overlaps=20
#                   ,nudge_x = lbld[sig_bonf==1,ndg]
#                   ,direction="y"
#                   ) + 
  theme(legend.position="bottom") + 
#  theme(text=element_text(face="bold")) + 
  coord_fixed() 

pw3p[,islbl:=as.numeric(target %in% lbld[["target"]])]
spatialp <- 
ggplot(pw3p[order(islbl)], aes(x=fold_change_v2
                 ,y=`0.5quant`,
                 fill=spatial_re_sd_mean, shape=)) + 
  geom_point(data=pw3p[p.value_v2 < 0.05/.N]
             ,aes(shape=8),size=3,stroke=1.2) + 
  geom_point(aes(shape=21),size=3) + 
  geom_point(data=pw3p[sig_bonf==1]
             ,aes(shape=0),size=6,stroke=1.2) + 
  scale_shape_identity(guide=guide_legend(override.aes=list(size=6))
                       ,breaks=c(0,8)
                       ,name=""
                       ,labels=c("9 Significant genes\n(INLA spatial RE model)"
                                 ,"114 Significant genes\n(non-spatial semi-naive model)")) + 
  theme_bw(base_size = 20) + 
  labs(x="Fold change (without spatial RE)"
       ,y ="Fold change (INLA, with spatial RE)") + 
#  scale_color_gradient2() + 
  scale_fill_gradient2(name = "Std Dev of\nSpatial\nRandom Effect"
                       ,high=scales::muted("red"), mid="white", low=scales::muted("blue")
                       #,midpoint = median(pw3p[["spatial_re_sd_mean"]])
                       ,midpoint = mean(pw3p[["spatial_re_sd_mean"]])
                       ) + 
  geom_abline(slope=1,intercept=0,lty=2) + 
  scale_x_continuous( n.breaks=20) + 
  scale_y_continuous( n.breaks=20) + 
  #geom_label_repel(data=lbld[sig_bonf!=1]
  #                 ,aes(label=target)
  #                 ,max.overlaps=20
  #                 ,nudge_x = lbld[sig_bonf!=1,ndg]
  #                 ) + 
  geom_label_repel(data=lbld
                   ,aes(label=target)
                   ,max.overlaps=20
                   ,nudge_x = lbld[,ndg]
                   ,direction="y"
                   ) + 
#  geom_label_repel(data=lbld[sig_bonf==1]
#                   ,aes(label=target)
#                   ,max.overlaps=20
#                   ,nudge_x = lbld[sig_bonf==1,ndg]
#                   ,direction="y"
#                   ) + 
  theme(legend.position="bottom") + 
  theme(axis.text.x=element_text(size=18)) + 
#  theme(text=element_text(face="bold")) + 
  coord_fixed() 

spatialp <- spatialp + 
  geom_hline(yintercept = 1,lty=2) + 
  geom_vline(xintercept = 1,lty=2)

pw1p[p.value < 0.05/.N]
pw1p[p.value < 0.05/.N][,.N,by=.(fold_change < 1)]
pw1p[p.value < 0.05/.N][,.N,by=.(fold_change < 1)][,prop:=N/sum(N)][]

pw2p[p.value_v2 < 0.05/.N]
pw2p[p.value_v2 < 0.05/.N][,.N,by=.(fold_change_v2 < 1)]
pw2p[p.value_v2 < 0.05/.N][,.N,by=.(fold_change_v2 < 1)][,prop:=N/sum(N)][]
pw2p[p.value_v2 < 0.05/.N][p.value_v1 > 0.05/.N]

#process <- fread("~/data/dan/decoding/wtx/Run1423_WTx12_Embryo/S2/fovprocessing.log")
#process[,idx:=1:.N,by=.(V3)]
#optim <- process[,elapsed:=V1 - shift(V1,1),by=.(V3)][!is.na(elapsed)]
#optim[,as.numeric(elapsed)/60]
#
#
#process2 <- fread("~/data/dan/decoding/wtx/Run1423_WTx12_Embryo/S2/fovprocessing_currdev.log")
#process2[,idx:=1:.N,by=.(V3)]
#currdev <- process2[,elapsed:=V1 - shift(V1,1),by=.(V3)][!is.na(elapsed)]
#currdev[,as.numeric(elapsed)/60]
#
#
#comp <- 
#merge(optim[,.(fov=V3, elapsed)]
#      ,currdev[,.(fov=V3, elapsed)]
#      ,by="fov", suffixes=c("_optim", "_current"))
#
#comp[,elapsed_diff:=elapsed_optim - elapsed_current]
#comp <- comp[,lapply(.SD, function(x) as.numeric(x) / 60),by=fov]
#comp[order(elapsed_current)]
#ggplot(comp, aes(x=elapsed_current, y=elapsed_optim)) + 
#  geom_point() + 
#  theme_bw() + 
#  geom_abline(slope=1,intercept=0,lty=2) + 
#  geom_smooth(method="lm", se=FALSE) + 
#  coord_fixed() + 
#  labs(x="elapsed optimized (minutes)", y="elapsed current mtm (minutes)", title="Elapsed mtm time per fov by version")
#

#> ms2[,head(.SD,1),by=.(target)][,summary(elapsed)]
#Min. 1st Qu.  Median    Mean 3rd Qu.    Max. 
#0.0630  0.1050  0.1240  0.1272  0.1490  0.3820 
#> ms1[,head(.SD,1),by=.(target)][,summary(elapsed)]
#Min. 1st Qu.  Median    Mean 3rd Qu.    Max. 
#0.0570  0.0810  0.1120  0.1124  0.1330  0.4470 
#> ms3[,head(.SD,1),by=.(target)][,summary(elapsed)]
#Min. 1st Qu.  Median    Mean 3rd Qu.    Max. 
#11.64   17.61   20.17   21.58   23.11  102.55 

#spatialp <- 
#ggplot(pw3p, aes(x=fold_change_v2
#                 ,y=`0.5quant`,
#                 fill=spatial_re_sd_mean)) + 
#  geom_point(data=pw3p[sig_bonf==1]
#             ,pch=0,size=6,stroke=1.2) + 
#  geom_point(data=pw3p[p.value_v2 < 0.05/.N]
#             ,pch=8,size=3,stroke=1.2) + 
#  geom_point(pch=21,size=3) + 
#  theme_bw(base_size=20) + 
#  labs(x="Fold change (without spatial RE)"
#       ,y = "Fold change (INLA, with spatial RE)") + 
##  scale_color_gradient2() + 
#  scale_fill_gradient2(name = "Std Dev of Spatial Random Effect"
#                       ,high=muted("red"), mid="white", low=muted("blue"), midpoint = median(pw2p[["contam_ratio"]])) + 
#  geom_abline(slope=1,intercept=0,lty=2) + 
#  scale_x_continuous( n.breaks=20) + 
#  scale_y_continuous( n.breaks=20) + 
#  geom_label_repel(data=lbld
#                   ,aes(label=target)
#                   ,max.overlaps=20
#                   ,nudge_x = lbld[,ndg]
#                   ) + 
#  theme(legend.position="bottom") + 
#  theme(text=element_text(face="bold")) + 
#  coord_fixed() 
#
#spatialp
#spatialp + 
#  scale_shape_manual(values = c(0, 8)
#                     ,labels=c("Significant with Spatial RE", "Significant w/out Spatial RE"))
#
pd <- merge(metainfo[,.(x_slide_mm, y_slide_mm,tissgrp)]
            ,sre[target == "WNT7A"]
            ,by.x=c("x_slide_mm", "y_slide_mm")
            ,by.y=c("x", "y"))#[order(mean)]
#pd <- copy(sre[target=="WNT7A"])
#setnames(pd, c("x", "y"), c("x_slide_mm", "y_slide_mm"))
wnt7ap <- 
ggplot(pd[order(mean)], aes(x_slide_mm, y_slide_mm, color=mean)) + 
  geom_point(pch=19,size=2) + 
  theme_bw(base_size=20) + 
  scale_color_gradient2(name = "Predicted\nSpatial\nRandom Effect"
                       ,high=scales::muted("red")
                       ,mid="white"
                       ,low=scales::muted("blue")
                       ,midpoint = mean(pd[["mean"]]))  + 
  coord_fixed()  + 
  theme(legend.position="bottom") + 
  labs(title="WNT7A Predicted\nSpatial Random Effects")

wnt7ap
pd <- merge(metainfo[,.(x_slide_mm, y_slide_mm,tissgrp)]
            ,sre[target == "CD79A"]
            ,by.x=c("x_slide_mm", "y_slide_mm")
            ,by.y=c("x", "y"))#[order(mean)]

cd79p <- 
ggplot(pd[order(mean)], aes(x_slide_mm, y_slide_mm, color=mean)) + 
  geom_point(pch=19,size=2) + 
  theme_bw(base_size=20) + 
  scale_color_gradient2(name = "Predicted\nSpatial\nRandom Effect"
                       ,high=scales::muted("red")
                       ,mid="white"
                       ,low=scales::muted("blue")
                       ,midpoint = median(pw2p[["contam_ratio"]]))  + 
  coord_fixed()  + 
  theme(legend.position="bottom", legend.text=element_text(size=8)) + 
  labs(title="CD79A Predicted\nSpatial Random Effects")


cowplot::plot_grid(spatialp + theme(plot.margin = unit(c(0.1,0,0,0.1), "cm")
                                    ,legend.box = "vertical"
                                   # ,legend.spacing.y = unit(0.0, "cm")
#                                   ,legend.key.size = unit(2, "cm")
                                    ) + 
#                     guides(fill=guide_colourbar(barwidth=10, "cm")) + 
                     labs(title="Impact of Spatial Random Effect Variation Fold Change Estimation and Inference ")
                   ,cd79p + theme(plot.margin = unit(c(0,0,0,0), "cm"))
                   ,wnt7ap + theme(plot.margin = unit(c(0,0.1,0,0), "cm"))
                   ,rel_widths = c(1, 0.5, 0.5)
                   ,labels = c("a", "b", "c")
                   ,nrow = 1) + 
  theme(plot.background = element_rect(fill = 'white',color='white')) + 
  cowplot::panel_border(remove=TRUE)


ggsave("colon_data/deplot3.png"
       ,width=1858*sqrt(25)
       ,height=919*sqrt(25)
       ,dpi=450
       ,units = "px")

plot3efg <- 
cowplot::plot_grid(spatialp + theme(plot.margin = unit(c(0.1,0,0,0.1), "cm")
                                    ,legend.box = "vertical"
                                    ,axis.text.x = element_text(size=12)
                                   # ,legend.spacing.y = unit(0.0, "cm")
#                                   ,legend.key.size = unit(2, "cm")
                                    ) + 
#                     guides(fill=guide_colourbar(barwidth=10, "cm")) + 
                     labs(title="Impact of Spatial Random Effect Variation\non Fold Change Estimation and Inference ")
                   ,cd79p + theme(plot.margin = unit(c(0.1,0,0,0), "cm")) 
                   ,wnt7ap + theme(plot.margin = unit(c(0.1,0.1,0,0), "cm"))
                   ,rel_widths = c(1, 0.5, 0.5)
                   ,labels = c("e", "f", "g")
                   ,label_size=20
                   ,nrow = 1) + 
  theme(plot.background = element_rect(fill = 'white',color='white')) + 
  cowplot::panel_border(remove=TRUE)

plot3efg






plist2 <- 
lapply(split(metainfo, by="tissgrp"), function(xx){
  ggplot(xx[clust!="B-cell"][order(ct)], aes(y_slide_mm, x_slide_mm, color=ct)) + 
    geom_point(size=0.05, alpha=0.5) + 
    geom_point(data=xx[clust=="B-cell"], size=0.3) + 
    coord_fixed() + 
    theme_bw() + 
    #scale_color_manual(values=rep(unname(pals::alphabet()),2)
    scale_color_manual(values=rev(ctpal2)
                       ,breaks=rev(names(ctpal2))
                       ,name="cell type"
                       ,guide=guide_legend(override.aes=list(size=4, alpha=1))) + 
    theme(legend.position="bottom")
    
})

leg1 <- cowplot::get_plot_component(plist[[1]] + 
                                      #theme(legend.position = "right")
                                      theme(legend.position = "right", legend.text = element_text(size=20), legend.title = element_text(size=20))
                                    , 'guide-box-right', return_all=TRUE)
upside <- 
cowplot::plot_grid(cowplot::plot_grid(plist2$ll + theme(legend.position="none") 
                                      ,cowplot::plot_grid(leg1
                                                          ,plist$ur + theme(legend.position="none")
                                                          ,nrow=2
                                                          ,rel_heights = c(1,3))
                                      ,nrow=1, rel_widths = c(6, 2.8))
                   ,labels=c("a", "")
                   ,label_size = 20) + 
  theme(plot.background = element_rect(fill = 'white',color='white')) + 
  cowplot::panel_border(remove=TRUE)
upside
ggsave("colon_data/colon_deplot_pt1.png"
       ,height =576*sqrt(30)
       ,width=1163*sqrt(30)
       ,units="px"
       ,dpi=450
       )
cowplot::plot_grid(upside, plot3cde, nrow=2, rel_heights = c(1, 1))




png("colon_data/deplot_3row.png"
    ,res=450
    ,width=19*sqrt(15)
    ,height=8*2.54*sqrt(15)
    ,units = "cm"
    )
cowplot::plot_grid(upside, pmiddle, plot3efg,nrow=3)
dev.off()
cowplot::plot_grid(upside, pmiddle,nrow=2)

png("colon_data/deplot_3row_v2.png"
    ,res=450
    ,width=19*sqrt(9)
    ,height=8*2.54*sqrt(9)
    ,units = "cm"
    )
cowplot::plot_grid(upside, pmiddle, plot3efg,nrow=3)
dev.off()

png("colon_data/deplot_3row_v3.png"
    ,res=450
    ,width=19*sqrt(4)
    ,height=8*2.54*sqrt(4)
    ,units = "cm"
    )
cowplot::plot_grid(upside, pmiddle, plot3efg,nrow=3)
dev.off()













cowplot::plot_grid(pmiddle,plot3efg,nrow=2)
ggsave("colon_data/deplot_2row.png"
       ,dpi=450
       ,width=19*sqrt(15)
       ,height=8*2.54/3*2 *sqrt(15)
       ,units = "cm"
       )




cowplot::plot_grid(upside, plot3cde, nrow=2, rel_heights = c(1, 1))
#ggsave("colon_data/deplottest.png"
#       ,height = 1482*20
#       ,width=1482*40
#       ,units="px"
#       ,dpi=450)

cowplot::plot_grid(leftside, rightside3, nrow=1)

leg2 <- cowplot::get_plot_component(plist_immcut[[1]] , 'guide-box-right', return_all=TRUE)


combleg <- cowplot::plot_grid(cowplot::ggdraw(leg1) + theme(plot.margin = unit(c(0,-5,0,-5), "cm"))
                              ,cowplot::ggdraw(leg2)+ theme(plot.margin = unit(c(0,-5,0,-5), "cm"))
                              , nrow=2)



ggsave("colon_data/deplot3cde.png"
       ,width=1858*sqrt(25)
       ,height=919*sqrt(25)
       ,dpi=450
       ,units = "px")


ggsave("colon_data/deplot3.png"
       ,width=1858*sqrt(25)
       ,height=919*sqrt(25)
       ,dpi=450
       ,units = "px")

hbbp <- 
ggplot(sre[target == "HBB"][order(mean)], aes(x, y, color=mean)) + 
  geom_point(pch=19,size=2) + 
  theme_bw() + 
  scale_color_gradient2(name = "Predicted Spatial Random Effect"
                       ,high=scales::muted("red")
                       ,mid="white"
                       ,low=scales::muted("blue")
                       ,midpoint = median(pw2p[["contam_ratio"]]))  + 
  coord_fixed()  + 
  theme(legend.position="bottom") + 
  labs(title="WNT7A Predicted Spatial Random Effects")

sre_sub <- merge(sre[target=="CCKBR"], metainfo_bcell, by.x=c("x", "y"), by.y=c("x_slide_mm", "y_slide_mm"))
ggplot(sre_sub[target == "CCKBR"][order(mean)], aes(x, y, color=mean)) + 
  geom_point(pch=19,size=2) + 
  theme_bw() + 
  scale_color_gradient2(name = "Predicted Spatial Random Effect"
                       ,high=scales::muted("red")
                       ,mid="white"
                       ,low=scales::muted("blue")
                       ,midpoint = median(pw2p[["contam_ratio"]]))  + 
#  coord_fixed()  + 
  theme(legend.position="bottom") + 
  labs(title="WNT7A Predicted Spatial Random Effects") + 
  facet_wrap(~tissgrp,scales="free")




ggplot(pw2p, aes(x=-log10(p.value_v2)
                 ,y=-log10(p.value_v1),
                 fill=contam_ratio)) + 
  geom_point(pch=21,size=3) + 
  theme_bw() + 
#  scale_color_gradient2() + 
  scale_fill_gradient2(high=muted("red"), mid="white", low=muted("blue"), midpoint = median(pw2p[["contam_ratio"]])
                       ,guide=guide_legend(override.aes=list(size=4))) + 
  geom_abline(slope=1,intercept=0,lty=2) + 
  scale_x_continuous(trans='log10', n.breaks=10) + 
  scale_y_continuous(trans='log10', n.breaks=10)  


smiDE::volcano(rbindlist(list(pw1[term=="immune_ct_neighbors"][,typ:="naive"]
                              ,pw2[term=="immune_ct_neighbors"][,typ:="v2"]
                              ))
               ,fold_change_cutoff = 2**0.5
               , interactive=FALSE
               ,limit_labels = TRUE)[[1]] + 
  facet_wrap(~typ,ncol=2)

pdcomb <- 
rbindlist(list(pw1[term=="immune_ct_neighbors"][,typ:="naive"]
                              ,pw2[term=="immune_ct_neighbors"][,typ:="v2"]
                              ))
smiDE::volcano(pdcomb
               ,fold_change_cutoff = 2**0.5
               , interactive=FALSE
               ,limit_labels = TRUE)[[1]] + 
  geom_line(data=pdcomb, aes(x=log2(fold_change), y= -log10(p.value), group=target)
            ,inherit.aes=FALSE)

ggplot(pw1, aes())
pw1[term=="immune_ct_neighbors"]


ggplot(pw1[term=="immune_ct_neighbors"], aes(x=log2(fold_change),y=-log10(p.value))) + 
  geom_point() + 
  scale_y_continuous(n.breaks=)


vnaive$`immune_ct_neighbors 14.54 / immune_ct_neighbors 9.458`

vp2 <- 
smiDE::volcano(pw2[term=="immune_ct_neighbors"]
               ,fold_change_cutoff = 2**0.5
               , interactive=FALSE
               ,limit_labels = TRUE)

vp2$`immune_ct_neighbors 14.54 / immune_ct_neighbors 9.458`

inlapw <- pw3[term=="immune_ct_neighbors"]
inlapw[,lcl_bonf:=`9.255831e-06quant`]
inlapw[,ucl_bonf:=`0.9999907quant`]
inlapw[,lcl:=`0.025quant`]
inlapw[,ucl:=`0.975quant`]

inlapw[,sig_bonf:=as.numeric(sign(log(lcl_bonf))==sign(log(ucl_bonf)))]
inlapw[,sig:=as.numeric(sign(log(lcl))==sign(log(ucl)))]
inlapw[,sig:=as.numeric(sign(log(`9.255831e-06quant`))==sign(log(`0.9999907quant`)))]


inlapw

inlapw[,.(sign(log(`9.255831e-06quant`)),sign(log(`0.9999907quant`)))]


vp3 <- 
smiDE::volcano(pw3[term=="immune_ct_neighbors"]
               ,fold_change_cutoff = 2**0.5
               , interactive=FALSE
               ,limit_labels = TRUE)
vp3$`immune_ct_neighbors 14.54 / immune_ct_neighbors 9.458`




pw
smiDE::volcano(pw[term=="immune_ct_neighbors"]
               ,fold_change_cutoff = 2^0.5
               ,interactive=FALSE
               ,limit_labels = TRUE
               )



dev3 <- readRDS("colon_data/de_results/de_v3.targset__0.rds")
smiDE::results(dev3, "pairwise")



##
ggplot(metainfo[clust %in% c("B-cell", other_immune_types)], aes(x_slide_mm, y_slide_mm, color=clust)) + 
  geom_point(size=0.1) + 
  theme_bw() + 
  coord_fixed() + 
  scale_color_manual(values=rep(unname(pals::alphabet2()),2) #
                     ,guide=guide_legend(override.aes=list(size=4))
  )

ggplot(metainfo_bcell[clust %in% c("B-cell", other_immune_types)], aes(x_slide_mm, y_slide_mm, color=clust)) + 
  geom_point(size=0.1) + 
  theme_bw() + 
  coord_fixed() + 
  scale_color_manual(values=rep(unname(pals::alphabet2()),2) #
                     ,guide=guide_legend(override.aes=list(size=4))
  )



celltype_nbrs <- 
  copy(pre_de_bcell$nblist$adjacency_counts_by_ct)[,other_immune_ct_neighbors:=rowSums(.SD)
                                                   ,.SDcols=c("B-cell", other_immune_types)
  ]



ggplot(celltype_nbrs[match(bcell_cellid, cell_ID)]
       ,aes(x=other_immune_ct_neighbors)) + 
  geom_histogram(bins = 21)

pd <- celltype_nbrs[match(bcell_cellid, cell_ID),hist(other_immune_ct_neighbors)]
pd <- celltype_nbrs[match(bcell_cellid, cell_ID),hist(other_immune_ct_neighbors)]

ggplot(pd, aes(x_slide_mm, y_slide_mm, color=other_immune_ct_neighbors)) + 
  geom_point(size=0.1) + 
  theme_bw() + 
  coord_fixed() 
#  scale_color_manual(values=rep(unname(pals::alphabet2()),2) #
#                     ,guide=guide_legend(override.aes=list(size=4))
#  )
