library(NeighVizGallery)

# install dependency
if(nchar(system.file(package='EBImage')) ==0){
  BiocManager::install("EBImage")
}


## define input and output folders
dataDir <- "File_path_to_folder_for_query_flow_cell"
main_out_dir <- "File_path_to_output_folder"
if(!file.exists(main_out_dir)) dir.create(main_out_dir)

## full cell metadata of the whole study
# one can get this from tertiary analysis 
# NOTE: cells not present in full cell metadata because of failed QC would still present in images
fullCellMeta <- read.csv(fs::path(main_out_dir, "allCells_metadata_coords.csv"))
rownames(fullCellMeta) <- fullCellMeta[['cell_ID']]


## example of choosing query cells ----
# pick query cells by random sampling 20 cells per cell type within given region 

# pick fov range 
fovLocs <- read.csv(fs::path(dataDir, "RunSummary", "latest.fovs.csv"), header = F)
colnames(fovLocs) <- c("slide", "x_stage_mm", "y_stage_mm", "z_stage_mm", "Zoffset", "ROI", "fov", "acqOrder")[1:ncol(fovLocs)]
fovToUse <- fovLocs[(fovLocs$x_stage_mm > -12.5) & (fovLocs$y_stage_mm  > 0), 'fov']

# pick cells per cell type
celltype_coln <- 'celltype'
cellNum_eachCT <- 20


tmpDF <- fullCellMeta[fullCellMeta$slide ==1 & fullCellMeta$fov %in% fovToUse, ]
tmpDF <- split(tmpDF, as.factor(tmpDF[[celltype_coln]]))

cells_to_use <- lapply(
  tmpDF, 
  function(df){
    # random sample 
    n <- min(nrow(df), cellNum_eachCT)
    
    return(df[sample(nrow(df), n), ])
  }
)
cells_to_use <- do.call(rbind, cells_to_use)

rownames(cells_to_use) <- cells_to_use[['cell_ID']]

# define how to group chosen cells in visualization 
groupVector <- setNames(cells_to_use[[celltype_coln]],
                        nm = cells_to_use[['cell_ID']])

## extract images from different kinds of inputs  ----
## Case 1: input is CellStatsDir with per FOV cell segmentation info ----
## extract cell segmentation files for each slide/flow cell
fileDF <- getCellSegFiles(CellStatsDir = fs::path(dataDir, "CellStatsDir"), 
                          # slideID should match what's used in fullCellMeta
                          slideID = 1)


# extract neighborhood of chosen cells from per FOV input 
imgArrayList <- extractNeighborhoodImgs(pixel_size_um = 0.120280945, 
                                        cropSize_um = 50,
                                        cells_to_use = cells_to_use$cell_ID, 
                                        morphFileDF = fileDF, 
                                        spatColns = c('CenterX', 'CenterY'), 
                                        # flag for do histograme equalization at FOV level for all channels in morphology image
                                        # set of FALSE if want to extract the original intensity value 
                                        equalizeFOV_flag = FALSE, 
                                        coreNum = 5)


# show cell_id around one query cell
colors_CellId <- displayCellId_singleImg(imgLabel = imgArrayList$labels[[1]], 
                                         query_cellID = names(imgArrayList$labels)[1],
                                         # seed for random color generation
                                         seed = 123)

## to print any display into file: 
# open a graphic device, display figure and then close
png(filename = fs::path(main_out_dir, "output_img.png"))
displayCellId_singleImg(imgLabel = imgArrayList$labels[[8]], 
                        query_cellID = names(imgArrayList$labels)[8],
                        # seed for random color generation
                        seed = 123)
dev.off()


# in case of protein dataset, further extract protein images
proteinFileDF <- getProteinDirFiles(ProteinDir = fs::path(dataDir, "ProteinDir"), 
                                    # slideID should match what's used in fullCellMeta
                                    slideID = 1)

proteinArrayList <- extractNeighborhoodImgs_fromProteinDir(pixel_size_um = 0.120280945, 
                                                           cropSize_um = 50,
                                                           cells_to_use = cells_to_use$cell_ID, 
                                                           # which targets and in what order to load
                                                           targets_to_load = c("VIM", "CD45", "E-cad", "Desmin"), 
                                                           proteinFileDF = proteinFileDF, 
                                                           spatColns = c('CenterX', 'CenterY'), 
                                                           # flag for do histograme equalization at FOV level for all channels in morphology image
                                                           # set of FALSE if want to extract the original intensity value 
                                                           equalizeFOV_flag = FALSE, 
                                                           coreNum = 5)
# the protein images could be combined with C902 morph images
allMorphArrayList <- combineMorphList(imgList1 = imgArrayList$morph,
                                      imgList2 = proteinArrayList$morph)


## Case 2: input is napari-CosMx stitched data and cellMeta contains fov and local coordiantes at per FOV level  ----
## this would extract stitched protein data too 

## extract cell segmentation files for each slide/flow cell from the stitched data 
napari_imgDir <- "File_path_to_images_folder_of_stitched_data"

stitchedInfo <- getGlobalCoords_forStitched(napari_imgDir, 
                                            cellMeta = fullCellMeta, 
                                            cellID_coln = 'cell_ID', 
                                            cellLabel_coln = 'CellId',
                                            fov_coln = 'fov',
                                            spatColns = c('CenterX', 'CenterY'))

# add label value and stitched coordinates into metadata
fullCellMeta <- cbind(fullCellMeta, 
                      stitchedInfo[['cellMeta']][fullCellMeta[['cell_ID']], c('labVal', 'global_x_px', 'global_y_px', 'global_x_mm', 'global_y_mm')])

# optional to swap the order or select only a few channels to stitch, must includes `labels` 
zarrDF <- stitchedInfo[['zarrDF']]
zarrDF <- zarrDF[c('labels', 'PanCK', 'Green', 'CD298_B2M', 'CD45', 'DNA'), ]

# extract neighborhood of chosen cells from stitched input
imgArrayList2 <- extractNeighborhoodImgs_fromStitched(
  pixel_size_um = stitchedInfo[['dimSetup']][['pixel_size_mm']]*1000, # pixel size in micron unit
  cropSize_um = 50,
  cells_to_use = cells_to_use[['cell_ID']], 
  zarrDF = zarrDF, 
  global_size = stitchedInfo[['dimSetup']][['global_size']], 
  cellMeta = fullCellMeta, 
  cellID_coln = 'cell_ID', 
  cellLabel_coln = 'labVal',
  spatColns = c('global_x_px', 'global_y_px'), 
  flip_axes = TRUE, # stitched images are in flipped orientation of per FOV images
  equalize_flag = FALSE,
  coreNum = 5)

# show cell_id around one query cell, with label values converted to cell_ID in legend 
colors_CellId2 <- displayCellId_singleImg(
  imgLabel = imgArrayList2$labels[[1]], 
  query_cellID = names(imgArrayList2$labels)[1],
  # seed for random color generation
  seed = 123,
  # show cell ID name instead of legend value 
  cellID_converter = setNames(fullCellMeta[['cell_ID']], 
                              nm = fullCellMeta[['labVal']])
)

## intensity adjustment of extracted image array ----
# One can increase the brightness of an image through addition, adjust the contrast through multiplication, and apply gamma correction through exponentiation.
# Let's change the intensity for 1st morph change of 5th cell
img1 <- imgArrayList2$morph[[5]][, , 1]
img1_comb <- EBImage::combine(
  img1, # original channel intensity 
  img1 + 0.3, # brighter
  img1 * 2, # higher contrast
  img1^0.5 # decrease gamma correction 
)
# show all 4 images in Plots window
EBImage::display(img1_comb, all = TRUE, method = 'raster')


## Different workflows of visualization and grouping list of image arrays into multi-frame cell gallery object ----
## Workflow option 1 ----
# create cell overlay on imageList first, then create multi-frame image object for each group, then display 

segList <- overlayCellBorder_imgList(rgbMorphImgs = imgArrayList$morph, 
                                     labelsImgs = imgArrayList$labels, 
                                     groupVector = groupVector, 
                                     morphToShow = c(4, 3, 1), # display 4th, 3rd, 1st channel in R, G, B color
                                     border_color = gplots::col2hex("white"), # border in white
                                     export_plots = 'none') # not to explot plots in png

# display 1st group in plots window
# show border on morph 
displayGallery_imgObj(imgObj = segList$combOverlays[[1]], 
                      cell_ids = segList$cell_IDs[[1]], 
                      imgType =  'morph', # overlay image is multi-channel
                      morphToShow = c(1,2,3), # show all 3 channels
                      method = 'raster', 
                      allFrame = TRUE,
                      colormode = 'color',
                      imagesPerRow = NULL)
# show cell labels only 
displayGallery_imgObj(imgObj = segList$combLabels[[1]], 
                      cell_ids = segList$cell_IDs[[1]], 
                      imgType =  'labels', # overlay image is multi-channel
                      morphToShow = c(1,2,3), # show all 3 channels
                      method = 'raster', 
                      allFrame = TRUE,
                      colormode = 'color',
                      imagesPerRow = NULL)

# display morph of 1st group in interactive viewer, show all channels in gray scale
displayGallery_imgObj(imgObj = segList$combRGBs[[1]], 
                      cell_ids = segList$cell_IDs[[1]], 
                      imgType =  'morph', # multi-channel
                      morphToShow = seq_len(dim(segList$combRGBs[[1]])[3]), # all 5 channels
                      method = 'browser', 
                      allFrame = TRUE,
                      colormode = 'grayscale',
                      imagesPerRow = NULL)


## Workflow option 2 ----
# create multi-frame image object for each group first, then get cell overlay and display

combLabels <- combineImgByGroup(imgList = imgArrayList$labels, 
                                groupVector = groupVector, 
                                imgType = 'labels')

combMorph <- combineImgByGroup(imgList = imgArrayList$morph, 
                               groupVector = groupVector, 
                               imgType = 'morph')

# overlay for 1st group and export all kind of plots as png 
segObj1 <- overlayCellBorder_imgObj(rgbMorphObj = combMorph$combImgObj[['astro.1']], 
                                    labelsObj = combLabels$combImgObj[['astro.1']], 
                                    cell_ids = combLabels$cell_IDs[['astro.1']], 
                                    morphToShow = c(4, 3, 1), # display 4th, 3rd, 1st channel in R, G, B color
                                    border_color =  '#ffffff', 
                                    export_plots = c('morph', 'labels', 'overlay'),
                                    plotresultdir = getwd(), 
                                    fileprefix = "group_astro1")

## create composite with other monochromes for morph
compObj1 <- createComposite_multiFrame(imgObj = combMorph$combImgObj[['astro.1']],
                                       vizChannels = c(4, 3, 1),  
                                       channelColors = c('Cyan', 'Magenta', 'Greys'), 
                                       equalizeFlags = NULL )  
# display new composite in plot window 
displayGallery_imgObj(imgObj = compObj1, 
                      cell_ids = combMorph$cell_IDs[['astro.1']], 
                      imgType =  'morph', # overlay image is multi-channel
                      morphToShow = c(1,2,3), # show all 3 channels
                      method = 'raster', 
                      allFrame = TRUE,
                      colormode = 'color',
                      imagesPerRow = NULL)



## color cells with metadata ----
# NOTE: cells without value in metadata would plot in same color as extracellular space in black 
# recommend to color cells before grouping, since the operation is done on per FOV level instead of per group level 

## Case 1: color by cell type
ctColors <- sort(unique(fullCellMeta[['refinedClusters_3.clust']]))
ctColors <- setNames(Seurat:::DiscretePalette(length(ctColors)+1, 
                                              palette = "polychrome")[-1], 
                     nm = ctColors)
ctMetaImgs <- colorCellsByMeta(labelsImgs = imgArrayList$labels, 
                               fullCellMeta = fullCellMeta, 
                               cellID_coln = 'cell_ID', 
                               cellLabel_coln = 'CellId',
                               # unique identifier for each fov, ok to use fov if only 1 slide
                               # if more than 1 slide, need to create additional column from slide_fov
                               UID_coln = 'fov', 
                               color_coln = 'refinedClusters_3.clust',
                               isDiscrete = TRUE,
                               # colors for each cell type; if NULL, determine automatically 
                               colors_to_use = ctColors, 
                               value_limits = NULL,
                               plotresultdir = getwd()) # to plot legend 

# group meta images as one multi-frame image object per group 
combMetaCTs <- combineImgByGroup(imgList = ctMetaImgs, 
                                 groupVector = groupVector, 
                                 imgType = 'morph') 
# display one group 
displayGallery_imgObj(imgObj = combMetaCTs$combImgObj[["L2"]], 
                      cell_ids = combMetaCTs$cell_IDs[["L2"]], 
                      imgType =  'morph', # metadata image is multi-channel
                      morphToShow = c(1,2,3), # show all 3 channels
                      method = 'raster', 
                      allFrame = TRUE,
                      colormode = 'color',
                      imagesPerRow = NULL)

# add cell borders to the meta images of one group 
segCTMeta1 <- overlayCellBorder_imgObj(rgbMorphObj = combMetaCTs$combImgObj[["L2"]], 
                                       labelsObj = combLabels$combImgObj[["L2"]], 
                                       cell_ids = combLabels$cell_IDs[["L2"]], 
                                       morphToShow = c(1,2,3), # display 3 channels in meta images
                                       border_color =  '#ffffff', 
                                       export_plots = 'none', # not to export in png
                                       plotresultdir = getwd(), 
                                       fileprefix = "ctMeta_L2")
# display the meta image with cell borders
displayGallery_imgObj(imgObj = segCTMeta1, 
                      cell_ids = combLabels$cell_IDs[["L2"]], 
                      imgType =  'morph', # metadata image is multi-channel
                      morphToShow = c(1,2,3), # show all 3 channels
                      method = 'raster', 
                      allFrame = TRUE,
                      colormode = 'color',
                      imagesPerRow = NULL)


## Case 2: color by numeric value, display 1~99% quantile range 
countMetaImgs <- colorCellsByMeta(labelsImgs = imgArrayList$labels, 
                                  fullCellMeta = fullCellMeta, 
                                  cellID_coln = 'cell_ID', 
                                  cellLabel_coln = 'CellId',
                                  # unique identifier for each fov, ok to use fov if only 1 slide
                                  # if more than 1 slide, need to create additional column from slide_fov
                                  UID_coln = 'fov', 
                                  color_coln = 'nCount_RNA',
                                  isDiscrete = FALSE, 
                                  # define seed colors for color bar; if NULL, use default Blue-Red divergent color scale
                                  colors_to_use = rev(RColorBrewer::brewer.pal(11, 'RdYlBu')), 
                                  # define the range of numeric values to display; if NULL, use 1~99% quantile range in entire study 
                                  value_limits = NULL,
                                  plotresultdir = getwd()) # to plot legend 


# group meta images as one multi-frame image object per group 
combMetaCount <- combineImgByGroup(imgList = countMetaImgs, 
                                   groupVector = groupVector, 
                                   imgType = 'morph') 
# display one group 
displayGallery_imgObj(imgObj = combMetaCount$combImgObj[["L2"]], 
                      cell_ids = combMetaCount$cell_IDs[["L2"]], 
                      imgType =  'morph', # metadata image is multi-channel
                      morphToShow = c(1,2,3), # show all 3 channels
                      method = 'raster', 
                      allFrame = TRUE,
                      colormode = 'color',
                      imagesPerRow = NULL)

# put metadata image on top of morphology image
countMetaMorph <- combineMetaToMorph_imgList(
  rgbMorphImgs = imgArrayList$morph, # a list of multi-channel morphology images 
  rgbMetaImgs = countMetaImgs,  # a list of 3-channel meta data images  
  morphToShow = c(4, 3, 1) # display 4th, 3rd, 1st channel of morphology in R, G, B color
)

# group and display one group
combMetaCountMorph <- combineImgByGroup(imgList = countMetaMorph, 
                                        groupVector = groupVector, 
                                        imgType = 'morph') 
displayGallery_imgObj(imgObj = combMetaCountMorph$combImgObj[["L2"]], 
                      cell_ids = combMetaCountMorph$cell_IDs[["L2"]], 
                      imgType =  'morph', # metadata image is multi-channel
                      morphToShow = c(1,2,3), # show all 3 channels
                      method = 'raster', 
                      allFrame = TRUE,
                      colormode = 'color',
                      imagesPerRow = NULL)

