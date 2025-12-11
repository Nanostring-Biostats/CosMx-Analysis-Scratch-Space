
### Function for 'total count/library size normalization' used for gaussian family ######
sparse_tc_norm2 <- function(sm){
  libsizes <- Matrix::colSums(sm)
  libsizes[libsizes==0] <- 1
  normed <- sm %*% Matrix::Diagonal(x=mean(libsizes)/libsizes)
  dimnames(normed) <- dimnames(sm)
  return(normed)
}

library(smiDE)
library(data.table)
RhpcBLASctl::omp_set_num_threads(1)  # limit high memory overuseage
RhpcBLASctl::blas_set_num_threads(1)

args <- commandArgs(trailingOnly=TRUE)
#args <- c("naiive", "ranef_tissue", "Treg", "poisson", "all")
#args <- c("contam_otherct", "ranef_tissue", "T_CD8_naive,T_CD8_memory", "gaussian", "all")

formtype <- args[1]
groupadj <- args[2]
ref_ct_arg <- args[3]
fam <- args[4]
slide_use <- args[5]
famchar <- fam


### Cell types runnign DE for.. combine T CD4 naive /memory into T CD4
###                             combine T CD8 naive/memory  into T CD8
### See example input argument above 
ref_ct <- gsub("T_C", "T C",ref_ct_arg)
ref_ct <- gsub("4_m", "4 m",ref_ct)
ref_ct <- gsub("4_n", "4 n",ref_ct)
ref_ct <- gsub("8_m", "8 m",ref_ct)
ref_ct <- gsub("8_n", "8 n",ref_ct)
ref_ct <- strsplit(ref_ct, ",")[[1]]


gem <- readRDS("SMI_Giotto_Object.rds")
metainfo <- copy(gem@cell_metadata$rna)

### Combine tumor clusters into a common cell type name 'tumor'
metainfo[,cell_type_all:=cell_type]
metainfo[grep("tumor", cell_type_all),cell_type_all:="tumor"]

### add spatial locations into the metadata
metainfo <- merge(metainfo, gem@spatial_locs$raw, by = "cell_ID",sort=FALSE)

formla <- switch(
  formtype 
  ,contam_otherct = ~niche + RankNorm(otherct_expr)
  ,naiive = ~niche
)

nbexpr_to_calc <- switch(
  formtype 
  ,contam_otherct = c("otherct")
  ,naiive = "otherct" 
)

normed <- gem@expression$rna$normalized
counts <- gem@expression$rna$raw
totalcounts <- metainfo$totalcounts
names(totalcounts) <- colnames(counts)
if(slide_use=="all"){
  metainfo <- metainfo[Run_Tissue_name %in% 
                         c("Lung5_Rep2", "Lung9_Rep2", "Lung6", "Lung13", "Lung12"
                           )] ## remove replicate slides, i.e., Lung5_Rep1, Lung5_Rep3, Lung9_Rep1
  formla <- switch(groupadj
                   ,ranef_tissue = update.formula(formla, ~.+ (1|Run_Tissue_name))
                   ,fixef_tissue = update.formula(formla, ~.+ Run_Tissue_name)
                   ,none = formla
  )
} else {
  metainfo <- metainfo[Run_Tissue_name %in% slide_use]
  normed <- normed[,metainfo$cell_ID]
  counts <- counts[,metainfo$cell_ID]
}

if(famchar!="gaussian") formla <- update.formula(formla, ~.+offset(log(totalcounts)))

##### 
system.time({
  pre_de_obj <- 
    smiDE::pre_de(metadata = metainfo
                  ,ref_celltype = ref_ct
                  ,cell_type_metadata_colname = "cell_type_all"
                  ,split_neighbors_by_colname = "Run_Tissue_name"
                  ,mm_radius = 0.05
                  ,sdimx_colname = "sdimx"
                  ,sdimy_colname = "sdimy"
                  ,cellid_colname = "cell_ID"
                  ,adjacencies_only = TRUE
    )
})
metainfo[cell_type %in% ref_ct,.N,by=Run_Tissue_name]

assay_mat <- counts 
if(famchar=="gaussian") assay_mat <- normed


de_obj <- 
  smiDE::smi_de(
    assay_matrix = assay_mat #assay_mat[,metainfo[cell_type_all %in% ref_ct][["cell_ID"]]]
    ,metadata = metainfo[cell_type_all %in% ref_ct]
    ,groupVar = "niche" 
    ,formula = formla
    ,family=famchar
    ,targets= rownames(counts)[1:5]#NULL #
    ,neighbor_expr_overlap_weight_colname = "weight"
    ,neighbor_expr_cell_type_metadata_colname = "cell_type_all"
    ,pre_de_obj = pre_de_obj
    ,neighbor_expr_overlap_agg = "sum"
    ,neighbor_expr_totalcount_scalefactor = totalcounts
    ,neighbor_expr_totalcount_normalize = ifelse(famchar=="gaussian", FALSE, TRUE)
    ,nCores=4
  )

saveRDS(de_obj, file = paste0("de_results/de"
                              ,".formtype__",formtype
                              ,".groupadj__",groupadj
                              ,".ref_ct__",ref_ct_arg
                              ,".fam__", fam
                              ,".slide_use__",slide_use
                              ,".rds")
)

