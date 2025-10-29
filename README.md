
# scPearsonPCA: An R package for efficient Principal Components decomposition of pearson residual normalized single cell data

## Introduction

Pearson residuals can be an effective normalization method for
spatially-resolved transcriptomic data. However, the normalized data
matrix produced is **dense** with nearly all values non-zero. This
limits the usability of pearson residuals for high-plex or large
datasets where computing memory can be a limiting factor.

However, the principal component (PC) decomposition of quasi-poisson
pearson residuals can be computed directly from a sparse counts matrix
and a few easy-to-calculate summary statistics. The R package
[scPearsonPCA](https://github.com/Nanostring-Biostats/CosMx-Analysis-Scratch-Space/tree/Main/_code/scPearsonPCA)
enables this PC decomposition, including the returning the PCs as if the
normalized data had been centered, scaled, and influential values
clipped according to common best practices for PC analysis.

In this post we demonstrate a few quick examples for the package, as
well as how to use the PCA results returned for common downstream
analyses like UMAP and unsupervised clustering.

## Installation

You can install `scPearsonPCA` using the remotes package

``` r
remotes::install_github("Nanostring-Biostats/CosMx-Analysis-Scratch-Space",
                         subdir = "_code/scPearsonPCA", ref = "Main")
```

## Most used functions

- `sparse_quasipoisson_pca_seurat()` This function takes a counts matrix
  as input and returns a PCA Seurat-style reduction by default. The PCA
  is based on quasi-poisson pearson residual normalization. Similar to a
  NormalizeData $\rightarrow$ ScaleData $\rightarrow$ RunPCA workflow,
  arguments available for returning the PCA as if data were centered and
  scaled, with extreme values clipped.

- `sparse_quasipoisson_pca_seurat_batch()` This function employs a
  batch-correction method to PCA, provided there is a single categorical
  variable for ‘batch’. Examples could be ‘patient’, ‘slide_id’,
  ‘tma_id’, etc. Akin to a quasipoisson glm-based correction with the
  batch variable used as a covariate, the basic idea here is that gene
  frequencies are computed separately by batch rather than across the
  full sample.

- `sparse_quasipoisson_pca_seurat_multiomic()` This function is useful
  for getting the PCA decomposition from same-cell multiomic data. Here
  we assume we have two matrices, a sparse RNA counts matrix `x`, and a
  (possibly dense,lower dimensional) matrix of already-normalized
  protein expression `xother`. We want to be able to get PCA from the
  combined dataset, while still applying effective normalization to the
  RNA data. This function uses similar algebraic shortcuts which allow
  for efficient calculation of PCA from the **covariance matrix** of
  pearson residuals from the RNA, by deriving PCA from the combined
  covariance matrix of pre-normalized protein, and pearson-normalized
  RNA. We’ll demo a full example below using a multiomics breast cancer
  dataset.

## Quickstart / Example Usage

First we’ll read in our Seurat object dataset. This is 960-plex NSCLC
tissue.

``` r
library(data.table); setDTthreads(1)
library(scPearsonPCA)
library(ggplot2)
library(ggrepel)
```

``` r
sem <- readRDS("seurat_obj_nsclc.rds")
```

## Pearson PCA using `sparse_quasipoisson_pca_seurat`

We typically recommend to use up to 2,000-3,000 highly variable genes
for generating PCA results.

A few summary stats used as input for calculating pearson-normalized
PCs:

- `tc` is the total counts per cell across all genes.
- `genefreq` is the proportion of transcript calls for each gene across
  all cells.

In general, it’s recommended to calculate these summary stats based on
all genes rather than a subset of highly variable genes. While this
dataset is only 960-plex, the code below uses 900 highly-variable genes
as a demonstration for a basic workflow.

``` r
tc <- Matrix::colSums(sem[["RNA"]]@counts) ## total counts per cell (across all genes)
genefreq <- scPearsonPCA::gene_frequency(sem[["RNA"]]@counts) ## gene frequency (across all cells)
sum(genefreq)==1 # TRUE

sem <- Seurat::FindVariableFeatures(sem, nfeatures = 900)
hvgs <- sem@assays$RNA@var.features

### Returns a Seurat-style DimReduc object with 
### hvgs x pcs feature loadings
### cells x pcs cell embeddings
### gene-length vector of the mean pearson of residuals
### gene-length vector of the standard deviation of pearson residuals
pcaobj <- 
sparse_quasipoisson_pca_seurat(sem[["RNA"]]@counts[hvgs,]
                               ,totalcounts = tc
                               ,grate = genefreq[hvgs]
                               ,scale.max = 10 ## PCs reflect clipping pearson residuals > 10 SDs above the mean pearson residual
                               ,do.scale = TRUE ## PCs reflect as if pearson residuals for each gene were scaled to have standard deviation=1
                               ,do.center = TRUE ## PCs reflect as if pearson residuals for each gene were centered to have mean=0
                               )
```

### Using the output for downstream analysis

Here we’ll take the `pcaobj` we created and add the “`DimReduc`” to our
seurat object. We can use this for downstream Seurat package analysis.
Then, we’ll create a UMAP from these PC results using a helper function
`make_umap`, which will also return a cells x cells nearest-neighbors
graph we’ll use for unsupervised clustering.

``` r
umapobj <- scPearsonPCA::make_umap(pcaobj)
sem[["pearsonpca"]] <- pcaobj$reduction.data
sem[["pearsonumap"]] <- umapobj$ump  ## umap
sem[["pearsongraph"]] <- Seurat::as.Graph(umapobj$grph) ## nearest neighbors / adjacency matrix used for unsupervised clustering
sem <- Seurat::FindClusters(sem, graph = "pearsongraph")
sem@meta.data$pearson_clusters <- sem@meta.data$seurat_clusters
umapplot <- scPearsonPCA::plot_umap(umapreduc = "pearsonumap", clustercol = "pearson_clusters", semuse = sem)
print(umapplot)
```

<img src="readmedata/umapplot.png" width="6000" />

## Batch corrected Pearson PCA using `sparse_quasipoisson_pca_seurat_batch`

The syntax is very similar if we want to correct for a batch variable,
in this case ‘patient’. A few extra arguments are the cell identifier
`cellid_colname`, batch identifier `batch_variable`, and a data.frame
`obs` which should include both of these columns.

We also re-compute gene frequency passing these same additional
arguments. This returns a genes x batches matrix of gene frequencies.

``` r
genefreq_batch <- 
gene_frequency(sem[["RNA"]]@counts, obs = data.table(sem@meta.data)[,.(cell_ID, patient)]
               ,cellid_colname = "cell_ID"
               ,batch_variable = "patient")

Matrix::colSums(genefreq_batch)

### Returns a Seurat-style DimReduc object with 
### hvgs x pcs feature loadings
### cells x pcs cell embeddings
### gene-length vector of the mean pearson of residuals
### gene-length vector of the standard deviation of pearson residuals
pcaobj_batch <- 
sparse_quasipoisson_pca_seurat_batch(sem[["RNA"]]@counts[hvgs,]
                                     ,totalcounts = tc
                                     ,grate = genefreq_batch[hvgs,] ##
                                     ,obs =data.table(sem@meta.data)[,.(cell_ID, patient)] 
                                     ,batch_variable = "patient"
                                     ,cellid_colname = "cell_ID"
                                     ,scale.max = 10
                                     ,do.scale = TRUE
                                     ,do.center = TRUE
                                     )
```

### Using the output for downstream analysis

Here again we’ll take the `pcaobj_batch` object and use it to create a
UMAP and perform unsupervised clustering.

``` r
umapobj_batch <- scPearsonPCA::make_umap(pcaobj_batch)
sem[["pearsonbatchpca"]] <- pcaobj_batch$reduction.data
sem[["pearsonbatchumap"]] <- umapobj_batch$ump
sem[["pearsonbatchgraph"]] <- Seurat::as.Graph(umapobj_batch$grph) ## nearest neighbors / adjacency matrix used for unsupervised clustering
sem <- Seurat::FindClusters(sem, graph = "pearsonbatchgraph")
sem@meta.data$pearson_clusters_batch <- sem@meta.data$seurat_clusters
umapplotbatch <- scPearsonPCA::plot_umap(umapreduc = "pearsonbatchumap", clustercol = "pearson_clusters_batch", semuse = sem)
print(umapplotbatch)
```

<img src="readmedata/umapplotbatch.png" width="6000" />

For this dataset, there are epithelial cell clusters which are fairly
specific to each patient. This is expected in cancer tissue and may not
warrant any batch correction.

For sake of demonstration here, however, we can visualize how the batch
variable ‘patient’ separates more on the UMAP before vs. after
batch-correction, where patient clusters have higher overlap in UMAP
space.

In case of known batch effects, multiple methods of correction may be
considered, and the best performing method may well vary from dataset to
dataset. One alternative recommended method which often works well is
‘Harmony’, described in this
[post](https://nanostring-biostats.github.io/CosMx-Analysis-Scratch-Space/posts/batchcorrection/).

``` r
umapplot_patient_batch <- scPearsonPCA::plot_umap(umapreduc = "pearsonbatchumap", clustercol = "patient", semuse = sem,alpha=0.1)
umapplot_patient <- scPearsonPCA::plot_umap(umapreduc = "pearsonumap", clustercol = "patient", semuse = sem,alpha=0.1)

cp <- 
cowplot::plot_grid(umapplot_patient_batch + labs(title = "UMAP using batch-corrected pearson PCs")
                   ,umapplot_patient + labs(title = "UMAP using pearson PCs")
                   ,nrow=1) + 
  theme(plot.background = element_rect(fill = 'white',color='white')) + 
  cowplot::panel_border(remove=TRUE)  
print(cp)
```

<img src="readmedata/umapplot_patient_compare.png" width="6000" />

## Pearson PCA for multiomic data using `sparse_quasipoisson_pca_seurat_multiomic`

Below we’ll use a breast cancer dataset with WTx RNA + 64-plex protein
to demo multiomic PCA and downstream analysis.

##### Reading in RNA and protein datasets

``` r
## Note, using a v3-style Seurat Object below, i.e., `sem[["RNA"]]@counts`
## For v5-style, the syntax would be like `sem[["RNA"]]$counts`
options(Seurat.object.assay.version = "v3")  

sem_protein <- readRDS("Breast/sem_protein.rds")
sem_rna <- readRDS("Breast/sem_rna.rds")
sem_rna <- subset(sem_rna, nCount_RNA >=100)
```

##### Normalizing the protein data

The `sparse_quasipoisson_pca_seurat_multiomic` function won’t do any
special normalization to the protein data, so the data should be
normalized going in. Here we’ll use pearson residuals from a gamma
regression model as a normalization method, although other reasonable
methods might be considered and left to the user’s disgression.

``` r
## normalize each protein taking pearson residuals from gamma generalized linear model.
pgamm <- scPearsonPCA::pearson_normalize_gamma_glm(as.matrix(sem_protein[["RNA"]]@counts)
                                                   ,upper_q_thresh = 0.99
                                                   ,lower_q_thresh = 0.01
                                                   )
### add a '_protein' suffix to the protein names - helps differentiate proteins and RNA's with the same name. 
rownames(pgamm) <- paste0(rownames(pgamm), "_protein")

pgamm[1:5,1:5] ## normalized protein data
```

                          c_1_100_1  c_1_100_10 c_1_100_100 c_1_100_1000
    4-1BB_protein        -0.1592977  0.03704303  0.48381341 -0.054332422
    B7-H3_protein        -0.2409877 -0.08788226  0.08696422  0.006506165
    Bcl-2_protein        -0.3037842 -0.26226294  0.26853206 -0.278620332
    Beta-catenin_protein -0.7617500 -0.06967635  0.14393940 -0.430273385
    CCR7_protein         -0.1174789 -0.21163185  0.24066117  5.283740671
                         c_1_100_1001
    4-1BB_protein          -0.0686042
    B7-H3_protein           0.8100425
    Bcl-2_protein          -0.2199310
    Beta-catenin_protein   -0.3508033
    CCR7_protein           -0.3366898

##### PCA analysis derived from the joint covariance matrix of the the proteins and RNAs

First we’ll identify a subset of highly variable genes from the RNA
data. Then , we’ll pass in the RNA counts matrix and normalized protein
matrix to get the PCA.

``` r
## select genes from RNA data to use for PCA analysis.
## Will take the top 2k features as well as any pre-defined immune celltyping markers from HieraType
sem_rna <- Seurat::FindVariableFeatures(sem_rna, nfeatures = 2000)
use_genes <- c(sem_rna@assays$RNA@var.features)


### total counts and gene frequency from RNA data computed across full dataset
tc <- Matrix::colSums(sem_rna[["RNA"]]@counts)
genefreq <- scPearsonPCA::gene_frequency(sem_rna[["RNA"]]@counts)


### Principal components derived from the combined covariance matrix of pearson residual normalized RNA and protein.
pcaobj_multiomics <- scPearsonPCA::sparse_quasipoisson_pca_seurat_multiomic(sem_rna[["RNA"]]@counts[use_genes,] ## raw counts matrix for hvgs
                                                                       ,totalcounts = tc
                                                                       ,grate = genefreq[use_genes]
                                                                       ,xother = pgamm ## normalized protein
                                                                       ,do.scale = TRUE
                                                                       ,do.center = TRUE
                                                                       ,scale.max = 10
                                                                       ,ncores = 1
)
```

    PC_ 1 
    Positive:  CD127_protein, pan-RAS_protein, p53_protein, EGFR_protein, TCF7_protein, NF-kB p65_protein, Channel-PanCK_protein, CTLA4_protein, EpCAM_protein, Beta-catenin_protein 
           LAMP1_protein, Ki-67_protein, AZGP1, FASN, SERHL2, DBI, KRT19, CD138_protein, PIP, NAMPT 
    Negative:  CD39_protein, CD16_protein, CD40_protein, CD45_protein, VIM, CD4_protein, SPARC, CD74, COL1A1, COL3A1 
           HLA-DRB1, COL1A2, STING_protein, HLA-DRA, CTSB, CST3, SMA_protein, APOE, Vimentin_protein, ICAM1_protein 
    PC_ 2 
    Positive:  KRT19, AFMID, SERHL2, KRT8, FOXP3_protein, p53_protein, KRT7, CYP1B1, IgD_protein, LAG3_protein 
           ALDOA, MLPH, XBP1, NDRG1, CD24, CD123_protein, PGK1, Her2_protein, Channel-PanCK_protein, CD56_protein 
    Negative:  H3C15, H4C14, H4C11, H3-7, H2BC4, H3C2, H4C12, H3C13, H2AC18, H2BC18 
           H3C10, H1-2, H3C7, H2BC9, H2AC16, H1-5, H2BC12L, H2AC11, H3C11, H3C8 
    PC_ 3 
    Positive:  CD34_protein, CD123_protein, CD31_protein, COL4A1, SPARC, IGFBP7, COL4A2, PLVAP, GZMB_protein, SMA_protein 
           CD56_protein, HSPG2, LAG3_protein, FOXP3_protein, COL3A1, COL1A2, COL1A1, IgD_protein, LAMA4, COL6A2 
    Negative:  HLA-DRB1, CD74, HLA-DRA, APOE, Channel-CD68_protein, CTSD, HLA-DR_protein, LYZ, CD11c_protein, HLA-DQB1 
           IFI30, CD68_protein, APOC1, HLA-DPA1, C1QC, GPNMB, CD45_protein, C1QB, HLA-DQA1, CTSB 
    PC_ 4 
    Positive:  PD-L2_protein, GITR_protein, FOXP3_protein, 4-1BB_protein, LAG3_protein, CD56_protein, GZMB_protein, COL1A1, iNOS_protein, IgD_protein 
           Her2_protein, COL1A2, COL3A1, AEBP1, PD-1_protein, GZMA_protein, LUM, CD11b_protein, FN1, COL12A1 
    Negative:  PLVAP, CD34_protein, CD31_protein, COL4A1, A2M, EGFL7, AQP1, FLT1, KDR, COL4A2 
           HSPG2, SHANK3, PODXL, PECAM1, ENG, CD34, MCAM, ECSCR, INSR, NHERF2 
    PC_ 5 
    Positive:  COL1A1, COL1A2, COL3A1, FN1, AEBP1, POSTN, LUM, COL6A3, COL12A1, C1S 
           THBS2, FBN1, SULF1, COL6A1, COL5A1, HTRA3, SFRP2, EpCAM_protein, Beta-catenin_protein, CTHRC1 
    Negative:  CD123_protein, GZMB_protein, FOXP3_protein, PD-L2_protein, LAG3_protein, iNOS_protein, CD56_protein, 4-1BB_protein, Her2_protein, GITR_protein 
           IgD_protein, CD31_protein, CD11b_protein, CD34_protein, CD68_protein, GZMA_protein, PD-1_protein, STING_protein, A2M, PLVAP 

### Using the output for downstream analysis

``` r
umapobj_multiomics <- scPearsonPCA::make_umap(pcaobj_multiomics)


sem_rna[["multiumap"]] <- umapobj_multiomics$ump
sem_rna[["multigrph"]] <- Seurat::as.Graph(umapobj_multiomics$grph[colnames(sem_rna),colnames(sem_rna)])
sem_rna <- Seurat::FindClusters(sem_rna, graph.name="multigrph")
sem_rna@meta.data$clusters_multi_unsup <- sem_rna@meta.data$seurat_clusters


xyp_multi_unsup <- xyplot("clusters_multi_unsup", metadata = sem_rna@meta.data, ptsize = 0.05) + coord_fixed()
print(xyp_multi_unsup)
```

![](readmedata/xyp_multi_unsup.png)

``` r
umapp <- plot_umap(umapreduc = "multiumap", clustercol = "clusters_multi_unsup", semuse = sem_rna)
print(umapp)
```

![](readmedata/umapp_multi.png)

#### RNA and protein marker heatmaps

Below we’ll look at a couple of marker heatmaps for these unsupervised
clusters. These functions below use the [HieraType
package](https://github.com/Nanostring-Biostats/CosMx-Analysis-Scratch-Space/tree/Main/_code/HieraType).

``` r
#remotes::install_github("Nanostring-Biostats/CosMx-Analysis-Scratch-Space", 
#                        subdir = "_code/HieraType", ref = "Main")
library(HieraType)

fctbl_rna <- clusterwise_foldchange_metrics(sem_rna[["RNA"]]@counts[use_genes,]
                                            ,totalcounts = tc
                                            ,metadata = sem_rna@meta.data
                                            ,cluster_column = "clusters_multi_unsup"
                                            )

hm_rna <- add_ncells(marker_heatmap(fctbl_rna)) + labs(title = "rna heatmap")
print(hm_rna)
```

<div class="wider-image">

<img src="./readmedata/hm_rna_multi.png" alt="Wide Image">

</div>

``` r
fctbl_protein <- HieraType::clusterwise_foldchange_metrics_protein(sem_protein[["RNA"]]@counts[,rownames(sem_rna@meta.data)]
                                                        ,metadata = sem_rna@meta.data
                                                        ,cluster_column = "clusters_multi_unsup"
                                                        ,propd = (pgamm[,rownames(sem_rna@meta.data)] > 1)
                                            )
hm_protein <- add_ncells(marker_heatmap(fctbl_protein)) + labs(title = "protein heatmap")
print(hm_protein)
```

<div class="wider-image">

<img src="./readmedata/hm_protein_multi.png" alt="Wide Image">

</div>
