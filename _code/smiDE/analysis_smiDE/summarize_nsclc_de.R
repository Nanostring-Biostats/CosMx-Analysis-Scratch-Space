
library(data.table); setDTthreads(2)
library(smiDE)

RhpcBLASctl::omp_set_num_threads(1)  # limit high memory overuseage
RhpcBLASctl::blas_set_num_threads(1)

gem <- readRDS("SMI_Giotto_Object.rds")
metainfo <- copy(gem@cell_metadata$rna)
metainfo[,cell_type_all:=cell_type]
metainfo[grep("tumor", cell_type_all),cell_type_all:="tumor"]
metainfo <- merge(metainfo, gem@spatial_locs$raw, by = "cell_ID",sort=FALSE)

##################################################
### Compute the overlap ratio metric (ORM)   #####
##################################################

pre_de_obj <- 
  smiDE::pre_de(counts = gem@expression$rna$raw
                ,normalized_data = gem@expression$rna$normalized
                ,metadata = metainfo
                ,ref_celltype ="all" 
                ,cell_type_metadata_colname = "cell_type_all"
                ,split_neighbors_by_colname = "Run_Tissue_name"
                ,mm_radius = 0.05
                ,sdimx_colname = "sdimx"
                ,sdimy_colname = "sdimy"
                ,weight_colname = "weight"
                ,aggregation = "sum"
                ,cellid_colname = "cell_ID"
                ,nb_celltypes_to_individually_calc = "otherct" 
                ,adjacencies_only = FALSE
  )

orm_metric <-
overlap_ratio_metric(assay_matrix = gem@expression$rna$normalized
                           ,adjacency_matrix = pre_de_obj$nblist$adjacency_mat
                           ,metada = metainfo
                           ,cluster_col = "cell_type_all"
                           )

orm_metric[ratio < 1,.N,by=.(cell_type_all)]
orm_metric[,cell_type:=cell_type_all]
orm_metric[grep("T CD8",cell_type_all),cell_type:="T_CD8_naive,T_CD8_memory"]
orm_metric[grep("T CD4",cell_type_all),cell_type:="T_CD4_memory,T_CD4_naive,Treg"]
orm_metric[,maxratio:=max(ratio),by=.(target,cell_type)]
fwrite(orm_metric, "nsclc_de_results/overlap_ratio_metric.csv")

###############################################################################################################
###############################################################################################################
######## Read in smiDE objects and combine results tables across DE scenarios          ##########
###############################################################################################################
###############################################################################################################


### combine the de outputs into .csv files for plotting and summarizing
def <- list.files("de_results_clean",pattern="de.*rds",full.names = TRUE)
ovrl <- pwl <- msl <- emml <- vector(mode='list',length=length(def))

for(ii in 1:length(def)){
  meth <- data.table(f=def[ii])
  meth[,frmla:=gsub(".groupadj__.*$", "", tstrsplit(f, "formtype__")[[2]])]
  meth[,groupadj:=gsub(".ref_ct__.*$", "", tstrsplit(f, "groupadj__")[[2]])]
  meth[,cell_type:=gsub(".fam__.*$", "", tstrsplit(f, "ref_ct__")[[2]])]
  meth[,fam:=gsub(".sand_se.*$", "", tstrsplit(f, "fam__")[[2]])]
  meth[,slide_use:=gsub(".rds", "", tstrsplit(f, "slide_use__")[[2]])]
  de <- readRDS(def[ii])
  ovrl[[ii]] <- cbind(meth, results(de, "one.vs.rest")[[1]]) ## one.vs.rest contrasts
  pwl[[ii]] <- cbind(meth, results(de, "pairwise")[[1]])[!grep("otherct_expr",term)] ## pairwise contrasts
  msl[[ii]] <- cbind(meth, results(de, "model_summary")[[1]]) ## model summaries
  emml[[ii]] <- cbind(meth, results(de, "emmeans")[[1]]) ## emmeans
  message(ii / length(def))
}
ovr <- rbindlist(ovrl,use.names=TRUE,fill=TRUE)
pw <- rbindlist(pwl,use.names=TRUE,fill=TRUE)
ms <- rbindlist(msl,use.names=TRUE,fill=TRUE)
emm <- rbindlist(emml,use.names=TRUE,fill=TRUE)


fwrite(ovr, file=paste0("plotdata/ovr.results.csv"))
fwrite(pw, file=paste0("plotdata/pw.results.csv"))
fwrite(ms, file=paste0("plotdata/ms.results.csv"))
fwrite(emm, file=paste0("plotdata/emm.results.csv"))

ovr[,f:=NULL]
pw[,f:=NULL]
ms[,f:=NULL]
emm[,f:=NULL]


################################################
### Attach overlap ratio output to DE outputs ##
################################################

## the column 'maxratio' is the column we'll use downstream

ovr <- fread("nsclc_de_results/ovr.results.csv") ## one.vs.rest contrast summaries
ms <- fread( "nsclc_de_results/ms.results.csv")  ## model summaries
emm <- fread("nsclc_de_results/emm.results.csv") ## emmeans summaries
ovr <- ovr[groupadj=="ranef_tissue"][sand_se=="FALSE"][slide_use=="all"][,-c("sand_se"),with=FALSE]
ms <- ms[groupadj=="ranef_tissue"][sand_se=="FALSE"][slide_use=="all"][,-c("sand_se"),with=FALSE]
#ovr[groupadj=="ranef_tissue"][sand_se=="FALSE"][slide_use=="all"][,.N,by=.(contrast,fam,cell_type,frmla)]
#ms[groupadj=="ranef_tissue"][slide_use=="all"][,.N,by=.(fam,cell_type,frmla,term)]

ovr <- 
merge(ovr, orm_metric[,head(.SD,1),by=.(target,cell_type)], by = c("target", "cell_type"),all.x=TRUE
      ,suffixes=c("", "_overlap"))

ms <- 
merge(ms, orm_metric[,head(.SD,1),by=.(target,cell_type)], by = c("target", "cell_type"),all.x=TRUE
      ,suffixes=c("", "_overlap"))


##################################################################################################
##################################################################################################
#                                                                                    #############
#  Summary stats to append:  get counts per cell by contrast level, slide, cell_type #############
#                                                                                    ############# 
##################################################################################################
##################################################################################################

include_slides <- c("Lung5_Rep2", "Lung9_Rep2", "Lung6", "Lung13", "Lung12")
cell_types <- ovr[,unique(cell_type)]
contrast_countsl <- 
  lapply(cell_types, function(xx){
    ct <- strsplit(xx, ",")[[1]]
    ct <- gsub("_", " ", ct)
    cells_count <- 
      metainfo[cell_type_all %in% ct][
        Run_Tissue_name %in% include_slides
      ]
      o <- data.table(stack(Matrix::rowSums(gem@expression$rna$raw[,cells_count[["cell_ID"]],drop=FALSE])))
      o <- o[,.(cell_type=xx,gene=ind, cell_type_cts = values, cell_type_n = nrow(cells_count))]
  })
contrast_counts_all <- rbindlist(contrast_countsl)
ovr <- merge(ovr
             ,contrast_counts_all[,.(cell_type,cell_type_n,cell_type_cts,target=gene)][,unique(.SD)]
             ,by=c("cell_type","target")
             ,all.x=TRUE,sort=FALSE)

ms <- merge(ms
             ,contrast_counts_all[,.(cell_type,cell_type_n,cell_type_cts,target=gene)][,unique(.SD)]
             ,by=c("cell_type","target")
             ,all.x=TRUE,sort=FALSE)


cell_type_cols_all <- c(RColorBrewer::brewer.pal(8,"Set1"))[1:7]
names(cell_type_cols_all) <- c("macrophage", "T CD4", "T CD8","Treg", "tumor","neutrophil", "plasmablast")
celltype_converter <- names(cell_type_cols_all)
names(celltype_converter) <- c("macrophage", "T_CD4_memory,T_CD4_naive,Treg", "T_CD8_naive,T_CD8_memory", "Treg", "tumor", "neutrophil", "plasmablast")
ovr[,cell_type_label:=celltype_converter[cell_type]]
ms[,cell_type_label:=celltype_converter[cell_type]]
ovr[,cpc_1:=counts_1 / ncells_1]
ovr[,cpc_cell_type:=cell_type_cts / cell_type_n]


###########################################################################
###                                                ########################
###  Attach Human protein atlas annotations        ########################
###                                                ########################
###########################################################################
imm <- fread("hpar_reference_data/hpar.immune_cells.csv")[variable=="Immune cell specificity"]
ctspec <- fread("hpar_reference_data/hpar.tissue_rna_expression.csv")[variable=="Single cell type specificity"]

ovr <- merge(ovr, imm[,.(target=symbol, imm_anno=value)],by="target",all.x=TRUE,sort=FALSE)
ms <- merge(ms, imm[,.(target=symbol, imm_anno=value)],by="target",all.x=TRUE,sort=FALSE)


#########################################################################
###                                                                    ##
###  Attach model messages / model-fitting time summaries from         ##
###                                                                    ##
#########################################################################

msgs <- ms[,.(target,cell_type,frmla,groupadj,fam,slide_use,model_msg=msg,elapsed)][,unique(.SD)]
ovr <- merge(ovr, msgs,all.x=TRUE,sort=FALSE)
ms <- merge(ms,   msgs,all.x=TRUE,sort=FALSE)

fwrite(ovr, file="nsclc_de_results/ovr.plotting.csv")
fwrite(ms, file="nsclc_de_results/ms.plotting.csv")

