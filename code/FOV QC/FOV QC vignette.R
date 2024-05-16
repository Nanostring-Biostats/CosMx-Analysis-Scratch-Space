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

# note: the correct value of max_prop_loss is a judgment call. Here we use a 
# small (conservative) value to produce the most instructive results; 
# for most purposes, we recommend using the default value of 0.3, 
# i.e. flagging FOVs where a bit has >30% signal loss.
res <- runFOVQC(counts = counts, xy = xy, fov = fov, barcodemap = barcodemap, max_prop_loss = 0.3) 
str(res)
# which FOVs were flagged:
res$flaggedfovs
# which genes are involved in the flagged bits in those FOVs:
head(res$flagged_fov_x_gene)
# list all the flagged genes (losing just one reporter cycle impacts a lot of genes):
res$flagged_fov_x_gene[, "gene"]

#### Explore results: -----------------

# Show which FOVs have been flagged:
mapFlaggedFOVs(res)

# Now let's dissect the causes of FOV flags. 
# First and more simply, look for loss in total signal strength:
FOVSignalLossSpatialPlot(res, shownames = TRUE) 
# Two FOVs with clear signal loss have been flagged; another with more marginal signal loss has passed. 

# heatmap of estimated bias suffered for each FOV * barcode bit (only flagged FOV * bits are colored);
# Use this to peek under the hood at the intermediate results used to flag individual FOVs.
FOVEffectsHeatmap(res) 
# We see FOV 19 has many flagged bits, though most are from just one color from a reporter cycle.
# Biological variability can cause single bits to be flagged, so we only flag FOVs where >=2
# However, all 4 bits of reporter cycle 12 are flagged in FOV 19, so we fail that FOV.
# Similarly, we fail FOV 17 for having 3 failed bits in both reporter cycle 18 and 27,
# and FOV 18 for failing all 4 bits of reporter cycle 18.
# FOVs 17 and 14 have sporadic bits flagged, but not consistently within a reporter cycle, so we pass them. 


## spatial plots of per-bit FOV effects:
# These plots show all bits from flagged reporter cycles (only reporter cycles 
# with >= 2 flagged bits in a single FOV get shown):
# Here we divide the tissue into sub-FOV grids, and we color each grid square by its change
# in a barcode bit's total gene expression from similar grid squares from other FOVs.
# FOVs flagged for the bit are highlighted yellow.
# Grid squares without enough information are omitted. 
par(mfrow = c(2,2))
FOVEffectsSpatialPlots(res = res, outdir = NULL, bits = "flagged_reportercycles") 
# Failures tend to be unambiguous: it's very clear when all 4 bits from a reporter
# cycle are lost in an impacted FOV.

# Below, see an example of a flagged bit that appears more driven by spatial biology than by FOV artifacts.
# No other bits from this reporter cycle have been flagged, so this result did not cause any FOVs to fail.
par(mfrow = c(1,1))
FOVEffectsSpatialPlots(res = res, outdir = NULL, bits = "reportercycle22R") 


