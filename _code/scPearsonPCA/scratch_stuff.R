make_umap <- function(pcaobj, min_dist=0.01, n_neighbors=30, metric="cosine",key ="UMAP_" ){
  ump <- 
    uwot::umap(pcaobj$reduction.data@cell.embeddings
               ,n_neighbors = n_neighbors
               ,nn_method = "annoy"
               ,metric = metric
               ,min_dist = min_dist
               ,ret_extra = c("fgraph","nn")
               ,verbose = TRUE)
  
  umpgraph <- ump$fgraph
  dimnames(umpgraph) <- list(rownames(ump$nn[[1]]$idx), rownames(ump$nn[[1]]$idx))
  #  umpgraph <- umpgraph + Matrix::Diagonal(n = length(rownames(ump$nn[[1]]$idx)), names = rownames(ump$nn[[1]]$idx))
  
  #colnames(ump$embedding) <- paste0(toupper(key), c(1,2)) 
  #Seurat::CreateDimReducObject(embeddings = ump$embedding, key = toupper(key))
  colnames(ump$embedding) <- paste0(key, c(1,2)) 
  ump <- Seurat::CreateDimReducObject(embeddings = ump$embedding, key = key)
  return(list(grph = umpgraph
              ,ump = ump))
}

sem[["pca"]] <- pcaobj$reduction.data
umapobj <- make_umap(pcaobj)

sem[["umap"]] <- umapobj$ump
sem[["umapgrph"]] <- Seurat::as.Graph(umapobj$grph)
sem <- Seurat::FindClusters(sem, graph.name = "umapgrph", resolution = 1.2)

md <- data.table(sem@meta.data)
md[["cd3d"]] <- sem[["RNA"]]@counts["CD3D",md[["cell_ID"]]]

qp_glm <- 
glm(cd3d ~ patient + offset(log(totalcounts)), data = md
    ,family = quasipoisson(link='log'))


summary(qp_glm)
exp(qp_glm$model)

qp_res <- residuals(qp_glm, type="pearson")
qp_res[1:10]


devtools::load_all("~/data/dan/scPearsonPCA")
pcabatch1 <- scPearsonPCA::sparse_quasipoisson_pca_seurat_batch(sem[["RNA"]]@counts
                                                               ,obs = sem@meta.data
                                                               ,batch_variable = "patient"
                                                               ,cellid_colname = "cell_ID"
                                                               ,scale.max = 10
                                                               ,do.scale = FALSE
                                                               ,do.center = FALSE
)


head()
scx <-  
  Matrix::t(ytilde) - 
  Matrix::Diagonal(x=sqrt(totalcounts)) %*% batch_mat  %*%
  Matrix::t(root_grate_phi_row)

plot(scx[,"CD3D"] , qp_res)
head(scx[,"CD3D"] / qp_res)
1/ sqrt(1.01)




abline(0,1)


source("~/data/dan/imputation_smoothing/markers_sets.R")
devtools::load_all("~/data/dan/imputation_smoothing/smiSmooth")
source("~/data/dan/smi_ecker/helper_functions.R")
sem <- readRDS("~/NAS_data/lwu/testRun/SpatialTest/giotto_test/SMITAP_proj/WTx_Exp32S1S2/AtoMx_flatFiles/seuratObject_BreastS1.RDS")

grate <- scPearsonPCA::gene_frequency(sem[["RNA"]]@counts)
sum(grate)
tc <- Matrix::colSums(sem[["RNA"]]@counts)
summary(tc)
sem@assays$RNA@var.features
sem <- Seurat::FindVariableFeatures(sem, nfeatures = 3000)
hvgs <- sem@assays$RNA@var.features
data.table(sem@meta.data)
useg <- c(hvgs, markers$immune, unlist(tcellmarkers))
useg <- intersect(useg, rownames(sem))

normed <- smiSmooth::totalcount_norm(Matrix::t(sem[["RNA"]]@counts[useg,])
                                     ,tc = tc)
fctbl_brst <- 
smiSmooth::clusterwise_foldchange_metrics(normed = Matrix::t(normed)
                   ,metadata = data.table(sem@meta.data)
                   ,cluster_column = "breastAtlas_semiClus_2")

setdiff(c("CD3D", "FOXP3", "CD8B", "CD4"), fctbl$gene)
setdiff(c("CD3D", "FOXP3", "CD8B", "CD4", "ITGAX", "LAMP3", "CLEC4C", "CLEC9A"), fctbl_brst$gene)
hm <- marker_heatmap(fctbl_brst, extras = c("CD3D", "FOXP3", "CD8B", "CD4", "ITGAX", "LAMP3"))
umapf(umapreduc = 'umap',semuse = sem, clustercol = "breastAtlas_semiClus_2")
fwrite(fctbl_brst, file="~/data/dan/imputation_smoothing/breast_reanalysis/fctbl_breast_semisup.csv")

pears <- 
scPearsonPCA::sparse_quasipoisson_pca_seurat(x = sem[["RNA"]]@counts[useg,]
                                             ,totalcounts = tc
                                             ,grate = grate[useg]
                                             ,scale.max = 10
                                             ,do.center = TRUE
                                             ,do.scale = TRUE
                                             ,ncores = 1)


saveRDS(pears, file="~/data/dan/imputation_smoothing/breast_reanalysis/pears_pca.rds")

sem[["pca"]] <- pears$reduction.data
umapobj <- make_umap(pears)

sem[["umap"]] <- umapobj$ump
sem[["umapgrph"]] <- Seurat::as.Graph(umapobj$grph)
sem <- Seurat::FindClusters(sem, graph.name = "umapgrph")
sem@meta.data$clusters_unsup <- sem@meta.data$seurat_clusters
umapf(umapreduc = 'umap',semuse = sem, clustercol = "clusters_unsup")


fc_unsup <- clusterwise_foldchange_metrics(normed = Matrix::t(normed), metadata = data.table(sem@meta.data), cluster_column = "clusters_unsup")
fwrite(fc_unsup, file="~/data/dan/imputation_smoothing/breast_reanalysis/fctbl_breast_unsup.csv")
marker_heatmap(fc_unsup, extras = c("CD3D", "FOXP3", "CD8B", "CD4", "ITGAX", "LAMP3"))


saveRDS(sem, file="~/data/dan/imputation_smoothing/breast_reanalysis/sem.addedto.rds")

sem@reductions
umap







