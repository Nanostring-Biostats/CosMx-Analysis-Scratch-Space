
library(data.table)
library(ggplot2)

load("colon_data/coloncancer.RData")

args <- commandArgs(trailingOnly = TRUE)

RhpcBLASctl::blas_set_num_threads(1)
metainfo <- data.table(cbind(annot, clust))

targset <- as.numeric(args[1])

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


bcell_cellid <- metainfo[clust=="B-cell", cell_ID]
pre_de_bcell$nblist$adjacency_counts_by_ct[match(bcell_cellid
                                                 ,cell_ID)
                                           ,summary(rowSums(.SD))
                                           ,.SDcols=setdiff(colnames(pre_de_bcell$nblist$adjacency_counts_by_ct), "cell_ID")]

immune_types <- c("B-cell", "macrophage", "mDC", "neutrophil"
                  ,"pDC", "mast", "monocyte", "NK", "plasmablast"
                  ,"T CD8 memory", "plasma.cell", "plasmablast"
                  ,"T CD4 naive", "Treg")

setdiff(res[,unique(first_type)], immune_types)

celltype_nbrs <- 
  copy(pre_de_bcell$nblist$adjacency_counts_by_ct)[,immune_ct_neighbors:=rowSums(.SD)
                                                   ,.SDcols=c(immune_types)
  ]





orm_metrics <- fread("colon_data/orm_metrics_colon.csv")
use_g <- orm_metrics[clust=="B-cell"][ratio < 1][["target"]]
targs <- data.table(g=use_g)[,tgrp:=floor((0:(.N-1)/200))]

normed <- Matrix::t(raw)
cs <- Matrix::colSums(normed)
normed <- normed %*% Matrix::Diagonal(x = mean(cs)/cs, names = colnames(normed))


metainfo_bcell[,`:=`(sdimx = x_slide_mm, sdimy = y_slide_mm)]
### pre de contam
pre_de_contam <- 
  smiDE::pre_de(counts = Matrix::t(raw)
                ,normalized_data = normed
                ,weight_colname = "weight"
                ,metadata = metainfo
                ,ref_celltype = "B-cell"
                ,cell_type_metadata_colname = "clust"
                ,sdimx_colname = "x_slide_mm"
                ,sdimy_colname = "y_slide_mm"
                ,mm_radius = 0.05 ### 15 microns
                ,contamination = "sum")


de_v3 <- 
  smiDE::smi_de(assay_matrix = Matrix::t(raw[bcell_cellid,])
                ,metadata = metainfo_bcell
                ,formula = ~immune_ct_neighbors  + RankNorm(otherct_expr) + offset(log(nCount_RNA))
                ,neighborhood_counts = pre_de_contam$nblist
                ,groupVar = "immune_ct_neighbors"
                ,family="nbinom2"
                ,spatial_model = list(name="GP_INLA"
                                      ,x_coord_col = "sdimx"
                                      ,y_coord_col = "sdimy"
                                      ,quantiles = c(0.025/nrow(targs), 0.025, 0.5, 0.975, 1-0.025/nrow(targs))
                                      ,num.threads=1
                                      
                )
                ,targets=targs[tgrp==targset][["g"]]
                ,nCores=2
  )


#results(de_v3, "pairwise", variable="immune_ct_neighbors")
#spatial_re_plot <- 
#  ggplot(results(de_v3, "spatial_random_effect")[[1]]
#         ,aes(x,y, fill=mean)) + 
#  geom_point(pch=21, color='black') + scale_fill_gradient2(midpoint = 0) + 
#  coord_fixed()  + 
#  facet_wrap(~target)


saveRDS(de_v3, file=paste0("colon_data/de_results/de"
                              ,".targset__", targset, ".rds")
        )

