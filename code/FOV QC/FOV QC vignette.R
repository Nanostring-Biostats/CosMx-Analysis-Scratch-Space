#### simple demo of how to use FOV QC:

## load barcodes file:
barcodemap <- readRDS("https://github.com/Nanostring-Biostats/CosMx-Analysis-Scratch-Space/code/FOV%20QC/6kbarcodes.RDS")
head(barcodemap)
str(barcodemap)

## load example data:
load("https://github.com/Nanostring-Biostats/CosMx-Analysis-Scratch-Space/code/FOV%20QC/FOV%20QC%20example%20data.RData")

# data structure:
str(counts)
counts[1:5, 1:5]
head(xy)
head(fov)

## run QC pipeline
res <- runFOVQC(counts = counts, xy = xy, fov = fov, barcodemap = barcodemap)
str(res)

## Explore results:

# heatmap of FOV * bit results
FOVEffectsHeatmap(res) 
# spatial plots of per-bit FOV effects:
dir.create("FOV_QC_per_bit_plots")
FOVEffectsSpatialPlots(res = res, outdir = "FOV_QC_per_bit_plots") 
# now you can examine the plots for various bits, especially those highlighted in the above heatmap. 
