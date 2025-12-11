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


targs <- data.table(g=colnames(raw))[,tgrp:=floor((0:(.N-1)/200))]


naive_de <- smiDE::smi_de(
  assay_matrix = Matrix::t(raw[bcell_cellid,])
  ,metadata = metainfo_bcell
  ,formula = ~immune_ct_neighbors + offset(log(nCount_RNA))
  ,groupVar = "immune_ct_neighbors"
  ,family="nbinom2"
  ,nCores = 2
  ,targets=targs[tgrp==targset][["g"]]
)

saveRDS(naive_de, file=paste0("colon_data/de_results/naive_de"
                              ,".targset__", targset, ".rds")
        )

