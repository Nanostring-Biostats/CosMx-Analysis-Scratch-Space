# CosMxDA Export Custom Module
message("Custom Script Version: 0.1")

# Copyright 2023 NanoString Technologies, Inc.
# This software and any associated files are distributed pursuant to the NanoString AtoMx Spatial 
# Information Platform Software as a Service Agreement, available on the NanoString Technologies, Inc 
# website at www.nanostring.com, as updated.  All rights reserved.  No permission is granted to modify, 
# publish, distribute, sublicense or sell copies of this software.

# Run the InSituCor package to explore spatial correlations in CosMx data.
# Return: 1. a .csv manifest of correlated gene modules
#         2. a RDS file of gene-gene correlations 
#         3. spatial plots of correlated modules

# User defined variables
# celltype      - Metadata column giving cell type. options: Leiden, Insitutype


#### advanced arguments not exposed to UI: --------------------------------------
return_all <- TRUE    # whether to return the full results of insitucor()
plot_res <- 600        # resolution for .png plots
verbose <- FALSE      # for debugging; prints progress if toggled
# For finer-grained control of insitucor(), the arguments can be changed in the function call.
# Just Control-F for "insitucor(" to find the relevant code. Run ?insitucor for guidance on the arguments. 

library(nanopipeline)
library(InSituCor)
library(RColorBrewer)
library(viridis)
library(pheatmap)
library(Matrix)
library(tiledb)
library(tiledbsc)


outDir <- "/output"
if (!dir.exists(outDir)) {
  dir.create(outDir)
}

if (verbose) {
  print("reached1")
  print(celltype)
  print(str(celltype))
}


#### Test User Variables ---------------------------------------------------
variableTest <- function(varName, varType, msg, required = TRUE){
  if(!varType %in% c("character", "logical", "numeric")){
    stop("varType must be \"character\", \"logical\", or \"numeric\"")
  }
  varMsg <- NULL
  typeClass <- switch(varType, 
                      character={"STRING"},
                      logical={"BOOL"},
                      numeric={"NUM"})
  if(required == TRUE){
    if(!exists(varName)){
      varMsg <- paste0("\n\"", varName, "\" ", "was not set in custom module creation")
    }
  }
  if(exists(varName)){
    if(varType == "logical"){
      if(is.null(get(varName))){
        assign(varName, FALSE, envir = .GlobalEnv)
      }
    }
    if(!is(get(varName), varType)){
      if(!is.null(get(varName)) | required == TRUE){
        varMsg <- paste0(varMsg, paste0("\n\"", varName, "\" varType was not set as ", typeClass, 
                                        " in custom module creation"))
      }
    }
    if(required == TRUE){
      if(is.null(get(varName))){
        varMsg <- paste0(varMsg, paste0("\n\"", varName, "\" was not set and is required"))
      }else{
        if(get(varName) == ""){
          varMsg <- paste0(varMsg, paste0("\n\"", varName, "\" was not set and is required")) 
        }
      }
    }
  }
  if(!is.null(varMsg)){
    varMsg <- paste0(msg, varMsg)
  }else{
    varMsg <- msg
  }
  return(varMsg)
}

variableMsg <- NULL
variableMsg <- variableTest(varName = "celltype", varType = "character", msg = variableMsg)

# stop on variable errors if applicable
if(!is.null(variableMsg)){
  stop(variableMsg)
}

if (verbose) {
  print("reached2")
}

#### extract basic data ---------------------------------------

# how to see column names: my_pipeline$tiledbsc_dataset$somas$RNA$obs$attrnames()

counts <- study$somas$RNA$X$members$counts$to_matrix()
counts <- Matrix::t(counts)
norm <- sweep(counts, 1, pmax(Matrix::rowSums(counts), 20), "/") * 500
rm(counts)

if (verbose) {
  print("reached2.4")
}

# align norm to metadata:
annotcellIDs <- study$somas$RNA$obs$to_dataframe("cell_ID")
annotcellIDs <- annotcellIDs$cell_ID

m <- match(annotcellIDs, rownames(norm))
norm <- norm[m, ]

# using slide ID as a stand-in for tissue, which won't generally be available
tissuevec <- study$somas$RNA$obs$to_dataframe("slide_ID_numeric")  
tissuevec <- tissuevec$slide_ID_numeric

# get xy:
xy <- as.matrix(cbind(study$somas$RNA$obs$to_dataframe("x_slide_mm"),
                      study$somas$RNA$obs$to_dataframe("y_slide_mm")))

if (verbose) {
  print("reached3")
  print(table(tissuevec))
  print(head(tissuevec))
  print(head(xy))
}



#### get cell type: ---------------------------------

# (dev note: this doesn't yet allow for the use of refineClusters, or for any
#  custom cell type assignments. This logic should eventually be expanded to do so.)

celltype <- match.arg(celltype, c("Insitutype", "Leiden"))
obs_colns <- study$somas$RNA$obs$attrnames()
if(celltype == "Insitutype"){
  celltypecolumn <- obs_colns[max(which(grepl("nbclust", obs_colns) & grepl("clusters", obs_colns)))] #(column on the right is the latest)
}
if (celltype == "Leiden") {
  celltypecolumn <- obs_colns[max(which(grepl("nn", obs_colns)))]
  #celltype <- grep("^nn_.*_cluster_cluster_", obs_colns, value = T)[1]
}
if(celltypecolumn < 0) {
  stop(paste0("couldn't find the cell type column for ", celltype))
}
celltypevec <- study$somas$RNA$obs$to_dataframe()[[celltypecolumn]]

if (verbose) {
  print("reached4")
  print(celltypecolumn)
  print(head(celltypevec))
  print(str(celltypevec))
  print(head(study$somas$RNA$obs$to_dataframe()))
  print(str(study$somas$RNA$obs$to_dataframe()))
  print(head(annotcellIDs))
  print(head(tissuevec))
  print(head(celltypevec))
  print(head(xy))
  print(norm[1:4,1:4])
}


# subset to complete data:
iscomplete <- as.logical((Matrix::rowSums(is.na(norm)) == 0) * !is.na(celltypevec) * !is.na(tissuevec) * (rowSums(is.na(xy)) == 0))
norm <- norm[iscomplete, ]
celltypevec <- celltypevec[iscomplete]
tissuevec <- tissuevec[iscomplete]
xy <- xy[iscomplete, ]


if (verbose) {
  print("reached5")
  print(table(iscomplete))
  print(dim(norm))
  print(length(celltypevec))
}

#### run InSituCor ---------------------------------------
# prep conditionon:
conditionon <-  data.frame(celltype = celltypevec)
rownames(conditionon) <- rownames(norm)

if (verbose) {
  print("reached6")
}

# remove genes with no information:
propnonzero <- colMeans(norm > 0)

# run insitucor:
res <- insitucor(       
  counts = norm[, propnonzero > 0.005], 
  conditionon = conditionon,
  celltype = celltypevec,
  neighbors = NULL,
  xy = xy,
  k = 50,
  radius = NULL,
  tissue = tissuevec,
  min_module_size = 3,
  max_module_size = 30,
  resolution = 0.02,
  corthresh = 0.1,
  min_module_cor = 0.1,
  gene_weighting_rule = "inverse_sqrt",
  roundcortozero = 0.01,
  max_cells = 5000,
  attribution_subset_size = 5000,
  verbose = TRUE
)

if (length(res) == 0) {
  stop("Insitucor failed to run")
}


if (verbose) {
  print("reached7")
}

#### save results ---------------------------------------

# return insitucor results list:
if (return_all) {
  # save smaller, diversely formatted results to download
  saveRDS(res$modules, file = paste0(outDir, "/InSituCor_modules.RDS"))
  saveRDS(res$condcor, file = paste0(outDir, "/InSituCor_conditional_correlation_matrix.RDS"))   
  saveRDS(res$attributionmats, file = paste0(outDir, "/InSituCor_celltype_involvement_per_module.RDS"))   
  saveRDS(res$celltypeinvolvement, file = paste0(outDir, "/InSituCor_celltype_involvement_all_modules.RDS"))   
} 


if (verbose) {
  print("reached8")
}

#### before plotting, arrange the tissues: -----------------------------


#' Arrange tissues in xy space to reduce whitespace
#' 
#' Uses a shelf algorithm: places tallest tissues on the bottom shelf, and so on. 
condenseTissues <- function(xy, tissue, tissueorder = NULL, buffer = 0.2, widthheightratio = 4/3) {
  
  # get each tissue's dimensions:
  tissdf <- data.frame(tissue = unique(tissue))
  tissdf$width <- sapply(unique(tissue), function(tiss) {
    diff(range(xy[tissue == tiss, 1], na.rm = T))
  })
  tissdf$height <- sapply(unique(tissue), function(tiss) {
    diff(range(xy[tissue == tiss, 2], na.rm = T))
  })
  
  # choose tissue order:
  if (!is.null(tissueorder)) {
    if (length(setdiff(tissdf$tissue, tissueorder)) > 0) {
      stop("values in tissue missing from tissueorder")
    }
    if (length(setdiff(tissueorder, tissdf$tissue)) > 0) {
      stop("values in tissueorder missing from tissue")
    }
    tissdf$order <- match(tissdf$tissue, tissueorder)
  } else {
    tissdf$order <- order(tissdf$height, decreasing = TRUE)
  }
  tissdf <- tissdf[tissdf$order, ]
  
  # choose number of tissues for first shelf:
  tissuesperrow <- min(round(sqrt(nrow(tissdf)) * widthheightratio * mean(tissdf$height) / mean(tissdf$width)), nrow(tissdf))
  targetwidth <- sum(tissdf$width[1:tissuesperrow]) + buffer * (tissuesperrow - 1)
  
  # place tissues:
  tissdf$x <- NA
  tissdf$y <- NA
  tempx <- 0
  tempy <- 0
  tempshelfheight <- 0
  tempshelfwidth <- 0
  
  for (i in 1:nrow(tissdf)) {
    
    # place this tissue:
    tissdf$x[i] <- tempx
    tissdf$y[i] <- tempy
    # update the shelf dimensions:
    tempshelfheight <- max(tempshelfheight, tissdf$height[i])
    tempshelfwidth <- tempx + tissdf$width[i]
    # move along the shelf:
    tempx <- tempx + tissdf$width[i] + buffer
    # start a new shelf if it's getting too wide:
    if (i < nrow(tissdf)) {
      if (abs(tempshelfwidth - targetwidth) < abs(tempshelfwidth + buffer + tissdf$width[i+1] - targetwidth)) {
        tempy <- tempy + tempshelfheight + buffer
        tempx <- 0
        tempshelfheight <- 0
        tempshelfwidth <- 0
      }
    }
  }
  # now update xy:
  for (tiss in unique(tissue)) {
    inds <- tissue == tiss
    xy[inds, 1] <- xy[inds, 1] - min(xy[inds, 1]) + tissdf$x[tissdf$tissue == tiss]
    xy[inds, 2] <- xy[inds, 2] - min(xy[inds, 2]) + tissdf$y[tissdf$tissue == tiss]
  }
  return(xy)  
}

if (verbose) {
  print("reached9.1")
  print(sum(is.na(xy)))
  print(sum(is.na(tissuevec)))
}

# first condense the xy space:
xy <- condenseTissues(xy = xy, 
                      tissue = tissuevec, 
                      tissueorder = NULL, 
                      buffer = 1, 
                      widthheightratio = 4/3) 

if (verbose) {
  print("reached9.2")
}

#### plotting ----------------------------------------------

#' Choose plot dimensions given a range in xy space:
#' Assumes points drawn with cex = 0.1
xy2inches <- function(span) {
  return(span * 0.7)  # well-considered scaling factor
}
# now choose plot size:
plotdims <- c(xy2inches(diff(range(xy[, 1]))), 
              xy2inches(diff(range(xy[, 2]))))

# plot module neighborhood scores:
sapply(colnames(res$scores_env), function(name) {
  modulegenes <- res$modules$gene[res$modules$module == name]
  png(paste0(outDir, "/InSituCor_module_", make.names(name), ".png"), width = plotdims[1], height = plotdims[2], res = plot_res, units = "in")
  par(mar = c(0,0,0,0))
  plot(xy, pch = 16, cex = 0.15, asp = 1,
       xaxt = "n", yaxt = "n", xlab = "", ylab = "",
       col = viridis_pal(option = "B")(101)[
         round(1 + 100 * pmin(res$scores_env[, name] / quantile(res$scores_env[, name], 0.99), 1))])
  for (tiss in unique(tissuevec)) {
    text(median(range(xy[tissuevec == tiss, 1], na.rm = T)),
         max(xy[tissuevec == tiss, 2], na.rm = T) + 0.2, 
         tiss, cex = 2)
    rect(min(xy[tissuevec == tiss, 1]), min(xy[tissuevec == tiss, 2]),
         max(xy[tissuevec == tiss, 1]), max(xy[tissuevec == tiss, 2]),
         lwd = 0.2, border = "black")
  }
  dev.off()
})

if (verbose) {
  print("reached10")
}

## plot attribution analysis results:
# modules x celltypes
pdf(paste0(outDir, "/InSituCor_attribution_moduleXcelltype.pdf"),
    width = pmin(2 + 0.2 * ncol(res$celltypeinvolvement), 12),
    height = pmin(2 + 0.2 * nrow(res$celltypeinvolvement), 20))
pheatmap(res$celltypeinvolvement, 
         col = colorRampPalette(c("white", "cornflowerblue"))(100), 
         main = "", fontsize_row = 8, fontsize_col = 8)
dev.off()

pdf(paste0(outDir, "/InSituCor_attribution_permodule_geneXcelltype.pdf"))
for (name in names(res$attributionmats)) {
  pheatmap(res$attributionmats[[name]], 
           col = colorRampPalette(c("white", "cornflowerblue"))(100),
           main = name, fontsize_row = 8, fontsize_col = 8)
}
dev.off()


if (verbose) {
  print("reached11")
}