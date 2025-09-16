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

setdiff(metainfo[,unique(clust)], immune_types)

celltype_nbrs <- 
  copy(pre_de_bcell$nblist$adjacency_counts_by_ct)[,immune_ct_neighbors:=rowSums(.SD)
                                                   ,.SDcols=c(immune_types)
  ]

metainfo_bcell <- merge(metainfo[clust=="B-cell"]
                        ,celltype_nbrs[,.(cell_ID, immune_ct_neighbors)]
                        , by="cell_ID")


orm_metrics <- fread("colon_de_results/orm_metrics_colon.csv")
use_g <- orm_metrics[clust=="B-cell"][ratio < 1][["target"]]
targs <- data.table(g=use_g)[,tgrp:=floor((0:(.N-1)/200))]

normed <- Matrix::t(raw)
cs <- Matrix::colSums(normed)
normed <- normed %*% Matrix::Diagonal(x = mean(cs)/cs, names = colnames(normed))

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
                ,aggregation = "sum"
                ,adjacencies_only = FALSE)

de_v2 <- smiDE::smi_de(
  assay_matrix = Matrix::t(raw[bcell_cellid,])
  ,metadata = metainfo_bcell
  ,formula = ~immune_ct_neighbors  + RankNorm(otherct_expr) + offset(log(nCount_RNA))
  ,groupVar = "immune_ct_neighbors"
  ,neighborhood_counts = pre_de_contam$nblist
  ,family="nbinom2"
  ,nCores = 2
  ,targets=targs[tgrp==targset][["g"]]
)

saveRDS(de_v2, file=paste0("colon_de_results/de_seminaive"
                              ,".targset__", targset, ".rds")
        )

