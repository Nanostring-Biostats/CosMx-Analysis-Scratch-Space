#### simple demo of how to use FOV QC:

## source the necessary functions:
source("https://raw.githubusercontent.com/Nanostring-Biostats/CosMx-Analysis-Scratch-Space/FOV-QC/code/FOV%20QC/FOV%20QC%20utils.R")

## load example data:
load(url("https://github.com/Nanostring-Biostats/CosMx-Analysis-Scratch-Space/raw/FOV-QC/code/FOV%20QC/FOV%20QC%20example%20data.RData"))

# data structure:
str(counts)
counts[1:5, 1:5]
head(xy)
head(fov)


#### run QC pipeline -------------------
res <- runFOVQC(counts = counts, xy = xy, fov = fov, barcodemap = barcodemap)
str(res)
# which FOVs were flagged:
res$flaggedfovs
# which genes are involved in the flagged bits in those FOVs:
head(res$flagged_fov_x_gene)

#### Explore results: -----------------

# heatmap of estimated bias suffered for each FOV * bit (flagged FOV * bits only):
FOVEffectsHeatmap(res) 

# spatial plots of per-bit FOV effects:
dir.create("FOV_QC_per_bit_plots")
FOVEffectsSpatialPlots(res = res, outdir = "FOV_QC_per_bit_plots", bits = "flagged") 
# now you can examine the plots for various bits, especially those highlighted in the above heatmap. 
