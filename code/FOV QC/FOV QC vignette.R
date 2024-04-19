#### simple demo of how to use FOV QC:

## source the necessary functions:
source("https://raw.githubusercontent.com/Nanostring-Biostats/CosMx-Analysis-Scratch-Space/FOV-QC/code/FOV%20QC/FOV%20QC%20utils.R")

## load barcodes:
allbarcodes <- readRDS(url("https://github.com/Nanostring-Biostats/CosMx-Analysis-Scratch-Space/raw/FOV-QC/code/FOV%20QC/barcodes_by_panel.RDS"))
names(allbarcodes)
# get the barcodes for the panel we want:
barcodemap <- allbarcodes$Hs_6k
head(barcodemap)

## load example data: a subset of FOVs and genes from a 6k panel study of breast cancer:
load(url("https://github.com/Nanostring-Biostats/CosMx-Analysis-Scratch-Space/raw/FOV-QC/code/FOV%20QC/FOV%20QC%20example%20data.RData"))

# data structure:
str(counts)
counts[1:5, 1:5]
head(xy)
head(fov)
# (Note : the above data objects can easily be extracted from a Giotto or Seurat or tileDB object. 
# See the docs for whichever system you're using for how to extract count matrices and metadata columns.
# Note that this code expects a count matrix with cells in rows & genes in columns, and 
# some toolkits store a transposed version of this matrix. In this case, use Matrix::t() to transpose
# your sparse counts matrix.)


#### run QC pipeline -------------------
res <- runFOVQC(counts = counts, xy = xy, fov = fov, barcodemap = barcodemap, max_prop_loss = 0.2)
str(res)
# which FOVs were flagged:
res$flaggedfovs
# which genes are involved in the flagged bits in those FOVs:
head(res$flagged_fov_x_gene)
# list all the flagged genes (losing just one reporter cycle impacts a lot of genes):
res$flagged_fov_x_gene[, "gene"]

#### Explore results: -----------------

# heatmap of estimated bias suffered for each FOV * barcode bit (only flagged FOV * bits are colored);
# Use this to peek under the hood at the intermediate results used to flag individual FOVs.
FOVEffectsHeatmap(res) 
# We see 2 bits (colors) flagged from 1 reporter cycle in 1 FOV, and 1 bit flagged in another FOV.
# Biological variability can cause single bits to be flagged, so we only flag FOVs where >=2
# bits from a single reporter cycle are flagged.

## spatial plots of per-bit FOV effects:
par(mar = c(1,1,3,1))
# show all bits from flagged reporter cycles (only reporter cycles with >= 2 flagged bits in a single FOV get shown):
# Here we divide the tissue into sub-FOV grids, and we color each grid square by its change
# in a barcode bit's total gene expression from similar grid squares from other FOVs.
# grid squares without enough information are omitted. 
par(mfrow = c(2,2))
FOVEffectsSpatialPlots(res = res, outdir = NULL, bits = "flagged_reportercycles") 
# it's clear that all 4 bits from this reporter cycle are lost in the impacted FOV.

# show all bits that were flagged:
par(mfrow = c(2,2))
FOVEffectsSpatialPlots(res = res, outdir = NULL, bits = "flagged_bits") 
# we can see that one bit (reportercycle18R) was flagged in a different FOV, but 
#  it appears to be biology-driven, not a FOV-level technical artifact
