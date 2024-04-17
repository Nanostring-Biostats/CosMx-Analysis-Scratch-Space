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


#### run QC pipeline -------------------
res <- runFOVQC(counts = counts, xy = xy, fov = fov, barcodemap = barcodemap, max_prop_loss = 0.2)
str(res)
# which FOVs were flagged:
res$flaggedfovs
# which genes are involved in the flagged bits in those FOVs:
head(res$flagged_fov_x_gene)

#### Explore results: -----------------

# heatmap of estimated bias suffered for each FOV * bit (only flagged FOV * bits colored):
FOVEffectsHeatmap(res) 

## spatial plots of per-bit FOV effects:
par(mar = c(1,1,3,1))
# show all bits from flagged reporter cycles (only reporter cycles with >= 2 flagged bits in a single FOV get shown):
par(mfrow = c(2,2))
FOVEffectsSpatialPlots(res = res, outdir = NULL, bits = "flagged_reportercycles") 

# show all bits that were flagged:
par(mfrow = c(2,2))
FOVEffectsSpatialPlots(res = res, outdir = NULL, bits = "flagged_bits") 
# we can see that one bit (reportercycle18R) was flagged in a different FOV, but it appears to be biology-driven, not artifactual
