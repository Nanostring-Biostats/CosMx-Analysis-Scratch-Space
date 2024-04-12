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

# heatmap of estimated bias suffered for each FOV * bit (flagged FOV * bits only):
FOVEffectsHeatmap(res) 

# spatial plots of per-bit FOV effects:
par(mar = c(1,1,3,1))
FOVEffectsSpatialPlots(res = res, outdir = NULL, bits = "flagged") 
# now you can examine the plots for various bits, especially those highlighted in the above heatmap. 
