#' @title getCellSegFiles
#' @description collect cell segmentation output file paths for CosMx assay 
#' @param CellStatsDir file path to CellStatsDir folder 
#' @param slideID slide ID that would be assigned to cells in the query flow cell
#' @return data.frame of file path for each fov's CellLabel, morphology 2D and cell stats files
#' @export 
getCellSegFiles <- function(CellStatsDir,
                            slideID = 1){
  CellStatsDir <- normalizePath(CellStatsDir)
  morphDir <- fs::path(CellStatsDir, "Morphology2D")
  
  # find morph images 
  if(!dir.exists(morphDir)){
    morphDir <- CellStatsDir
  }
  tmpFiles <- dir(morphDir, recursive = F, 
                  pattern = "S[0-9]_C902_P99_N99_F[0-9]+.TIF$", full.names = T)
  
  if(length(tmpFiles)<1){
    stop(message("No C902 TIF images are found in morphology folder: ", morphDir))
  }else{
    message(sprintf("Found %d C902 TIF images in folder: %s ", length(tmpFiles), morphDir))
  }
  
  fileDF <- lapply(
    tmpFiles, 
    function(x){
      data.frame(file_path = x, 
                 fov = as.numeric(gsub("F","", unlist(stringr::str_extract(basename(x), "F[0-9]+")))))
    }
  )
  fileDF <- do.call(rbind, fileDF)
  colnames(fileDF)[1] <- 'imgPath'
  if(sum(duplicated(fileDF[['fov']]))>0){
    stop("Multiple files found for same fov!")
  }
  
  morphFileDF <- fileDF 
  morphFileDF[['slide']] <- slideID
  
  
  
  # cell label images, search root folder if no morph2D folder
  if (morphDir == CellStatsDir){
    # old version of folder structure, all cell segmentation files in root folder 
    tmpFiles <- dir(CellStatsDir, recursive = F, 
                    pattern = "CellLabels_F[0-9]+.tif$", full.names = T)
  }else{
    # new version of folder structure, all cell segmenation files in FOV subfolder 
    tmpFiles <- dir(CellStatsDir, recursive = T, 
                    pattern = "CellLabels_F[0-9]+.tif$", full.names = T)
  }
  
  message(sprintf("Found %d CellLabels tif images.", length(tmpFiles)))
  
  fileDF <- lapply(
    tmpFiles, 
    function(x){
      data.frame(file_path = x, 
                 fov = as.numeric(gsub("F","", unlist(stringr::str_extract(basename(x), "F[0-9]+")))))
    }
  )
  fileDF <- do.call(rbind, fileDF)
  colnames(fileDF)[1] <- 'labelPath'
  if(sum(duplicated(fileDF[['fov']]))>0){
    stop("Multiple files found for same fov!")
  }
  morphFileDF <- merge(morphFileDF, fileDF)
  
  
  
  # cell stats file with local pixel centroid for each cell, search root folder if no morph2D folder
  if (morphDir == CellStatsDir){
    # old version of folder structure, all cell segmentation files in root folder 
    tmpFiles <- dir(CellStatsDir, recursive = F, 
                    pattern = "S[0-9]_Cell_Stats_F[0-9]+.csv$", full.names = T)
  }else{
    # new version of folder structure, all cell segmenation files in FOV subfolder 
    tmpFiles <- dir(CellStatsDir, recursive = T, 
                    pattern = "S[0-9]_Cell_Stats_F[0-9]+.csv$", full.names = T)
  }
  
  message(sprintf("Found %d Cell_Stats csv files.", length(tmpFiles)))
  
  fileDF <- lapply(
    tmpFiles, 
    function(x){
      data.frame(file_path = x, 
                 fov = as.numeric(gsub("F","", unlist(stringr::str_extract(basename(x), "F[0-9]+")))))
    }
  )
  fileDF <- do.call(rbind, fileDF)
  colnames(fileDF)[1] <- 'cellStatsPath'
  if(sum(duplicated(fileDF[['fov']]))>0){
    stop("Multiple files found for same fov!")
  }
  morphFileDF <- merge(morphFileDF, fileDF)
  
  message(sprintf("%d FOVs with all files for C902, CellLabels and Cell_Stats.", nrow(morphFileDF)))
  
  return(morphFileDF)
  
}



#' @title extractNeighborhoodImgs
#' @description extract the neighborhood images around query cells
#' @param pixel_size_um pixel size in um unit of raw images
#' @param cropSize_um dimension of extracted neighborhood for query cells in um unit
#' @param cells_to_use vector of cell_ID for query cells in `c_[slide]_[fov]_[CellId]` format 
#' @param morphFileDF data.frame of file path for the paired cell segmentation outputs of each fov
#' @param spatColns a vector of column names for xy spatial coordinates of each cell in Cell_Stats files, in pixel unit 
#' @param equalizeFOV_flag flag to equalize intensity of morph images for each channel at fov level (default = FALSE)
#' @param coreNum number of cores for mclapply parallelization 
#' @param label_bit encoded range for label image (default = 16)
#' @return a nested list of 2 elements, `labels` and `morph`, which contains list of cell-wise image array
#' @export 
extractNeighborhoodImgs <- function(pixel_size_um = 0.120280945, 
                                    cropSize_um = 50,
                                    cells_to_use, 
                                    morphFileDF, 
                                    spatColns = c('CenterX', 'CenterY'), 
                                    equalizeFOV_flag = FALSE, 
                                    coreNum = 5, 
                                    label_bit = c(16, 8, 32)){
  label_bit <- match.arg(as.character(label_bit), as.character(c(16, 8, 32)))
  label_bit <- as.numeric(label_bit)
  
  colns <- c('slide', 'fov', 'imgPath', 'labelPath', 'cellStatsPath')
  if(!all(colns %in% colnames(morphFileDF))){
    stop(sprintf("Missing essential columns in `morphFileDF`: `%s`.", 
                 paste0(setdiff(colns, colnames(morphFileDF)), collapse = "`, `")))
  }
  
  if(length(spatColns) !=2){
    stop("Must provide 2 spatail column names in `spatColns`.")
  }
  
  cropSize_px <- round(cropSize_um/pixel_size_um)
  message(sprintf("Extract %d px X %d px box around each query cells.", 
                  cropSize_px+1, cropSize_px+1))
  
  
  queryDF <- lapply(
    cells_to_use, 
    function(x){
      x2 <- as.numeric(strsplit(x, '_')[[1]])
      tt <- data.frame(cell_ID= x, 
                       slide = x2[2], 
                       fov = x2[3], 
                       CellId = x2[4])
      return(tt)
    })
  queryDF <- do.call(rbind, queryDF)
  queryDF <- merge(queryDF, morphFileDF[, colns, drop = F])
  
  if(nrow(queryDF)<1){
    stop("Found no matched record between `cells_to_use` and `morphFileDF`, please check if cell id is in `c_[slide]_[fov]_[CellId]` format.")
  }else{
    message(sprintf("Found record for %d query cells within %d fov.", nrow(queryDF), 
                    nrow(unique(queryDF[, c('slide', 'fov')]))))
  }
  
  queryDF[['UID']] <- paste0(queryDF[['slide']], '_', queryDF[['fov']])
  
  
  # get px coordinates for each cell
  tmpDF <- split(queryDF, as.factor(queryDF[['UID']]))
  tmpDF <- lapply(
    tmpDF, 
    function(df){
      cellStatsDF <- read.csv(df[1, 'cellStatsPath'])
      
      if(!all(c(spatColns, 'CellId') %in% colnames(cellStatsDF))){
        stop(sprintf("Missing essential columns (`%s`) in cell stats file: `%s`.", 
                     paste0(setdiff(c(spatColns, 'CellId'), colnames(cellStatsDF)), collapse = "`, `"),
                     df[1, 'cellStatsPath']))
      }
      
      cellStatsDF <- cellStatsDF[, c('CellId', spatColns)]
      cellStatsDF <- merge(df, cellStatsDF, by = 'CellId')
      
      return(cellStatsDF)
    }
  )
  tmpDF <- do.call(rbind, tmpDF)
  if(nrow(tmpDF) != nrow(queryDF)){
    message(spritf("%d out of %d query cells have NO cell stats information.", nrow(queryDF) - nrow(tmpDF), nrow(queryDF)))
  }
  
  queryDF <- tmpDF
  
  
  # inputs: df[, c('fov', 'imgPath','labelPath','cell_ID','CenterX','CenterY')], 
  # global values: cropSize_px, spatColns, equalizeFOV_flag
  # fill 0 for crop window smaller than others 
  extractPerFovFun <- function(df){
    uid <- df$uid[1]
    
    message('UID (slide_FOV) ', uid)
    
    # get morph image 
    img1 <- EBImage::readImage(df[1, 'imgPath'])
    fov_size <- dim(img1)[1:2]
    
    
    # get topLeft corner of crop box, clip at 1
    tmpCoord <- df[, spatColns] - floor(cropSize_px/2)
    colnames(tmpCoord) <- paste0('min_', c('x_px', 'y_px'))
    tmpCoord <- sweep(as.matrix(tmpCoord), 2, c(1,1), pmax)
    
    spatCoord <- cbind(df[, c('cell_ID', spatColns)], as.data.frame(tmpCoord))
    
    # get bottomRight corner, clip at maximum fov size 
    tmpCoord <- df[, spatColns] + ceiling(cropSize_px/2)
    colnames(tmpCoord) <- paste0('max_', c('x_px', 'y_px'))
    tmpCoord <- sweep(as.matrix(tmpCoord), 2, fov_size, pmin)
    
    spatCoord <- cbind(spatCoord, as.data.frame(tmpCoord))
    rownames(spatCoord) <- spatCoord[['cell_ID']]
    
    # equalize across channels for morph image 
    if(equalizeFOV_flag){
      img1 <- EBImage::equalize(img1, levels = 2**16)
    }
    
    # get label image
    img2 <- EBImage::readImage(df[1, 'labelPath'])
    # convert label to integer
    img2 <- round(img2*(2**label_bit))
    
    # combine channel together 
    img1 <- EBImage::abind(img1, img2, along = 3)
    
    # crop images for each cell
    img2 <- lapply(
      seq_len(nrow(spatCoord)), 
      function(idx){
        # return the image array by itself
        img1_crop <- img1[seq(spatCoord[idx, 'min_x_px'], 
                              spatCoord[idx, 'max_x_px']), 
                          seq(spatCoord[idx, 'min_y_px'], 
                              spatCoord[idx, 'max_y_px']), ]
        img1_crop <- img1_crop@.Data
        
        # fill 0 for region smaller than crop window 
        missX <- cropSize_px +1 - dim(img1_crop)[1]
        missY <- cropSize_px +1 - dim(img1_crop)[2]
        
        if(missX >0){
          mockData <- array(0, dim = c(missX, dim(img1_crop)[2:3]))
          img1_crop <- abind::abind(img1_crop, mockData, along = 1)
        }
        
        if(missY >0){
          mockData <- array(0, dim = c(dim(img1_crop)[1], missY, dim(img1_crop)[3]))
          img1_crop <- abind::abind(img1_crop, mockData, along = 2)
        }
        
        # split array into 2 elements 
        res <- list(labels = img1_crop[, , dim(img1)[3]], 
                    morph = img1_crop[, , seq_len(dim(img1)[3]-1)])
        
        return(res)
      })
    names(img2) <- rownames(spatCoord)
    
    return(img2)
  }
  
  # extract info for each fov 
  tmpDF <- split(queryDF, as.factor(queryDF[['UID']]))
  
  
  imgData <- parallel::mclapply(tmpDF, 
                                mc.cores = coreNum, 
                                FUN = extractPerFovFun)
  
  # get cell_id
  imgListCellIds <- unlist(lapply(imgData, function(x) names(x)))
  
  imgData <- do.call(c, imgData)
  names(imgData) <- imgListCellIds
  
  # restructure to return 2 list for each element separately 
  finalRes <- list(labels = lapply(imgData, '[[', 'labels'), 
                   morph = lapply(imgData, '[[', 'morph'))
  
  return(finalRes)
} 


#' @title getProteinDirFiles
#' @description collect protein image output file paths for CosMx assay 
#' @param ProteinDir file path to ProteinDir folder 
#' @param slideID slide ID that would be assigned to cells in the query flow cell
#' @return data.frame of file path for each fov's decoded protein images, and per cell stats files
#' @export 
getProteinDirFiles <- function(ProteinDir,
                               slideID = 1){
  ProteinDir <- normalizePath(ProteinDir)
  
  # fov folders
  fovDirs <- list.dirs(ProteinDir, full.names = F, recursive = F)
  fovDirs <- grep("FOV[0-9]+", fovDirs, value = T)
  fovDirs <- setNames(as.numeric(gsub("FOV", "", fovDirs)), 
                      nm = fovDirs)
  if(length(fovDirs)<1){
    stop(sprintf("Did not find any FOV folders under the provided `ProteinDir`: `%s`", 
                 ProteinDir))
  }else{
    message(sprintf("Found %d FOV folders.", length(fovDirs)))
  }
  
  
  # find ProteinImages and PerCellStats under each fovDirs
  allFileDF <- lapply(
    seq_along(fovDirs), 
    function(idx){
      # protein images
      tmpFiles <- dir(fs::path(ProteinDir, names(fovDirs)[idx], "ProteinImages"), 
                      pattern = "S[0-9]_C[0-9]+_F[0-9]+_.*.TIF$", recursive = F, full.names = T)
      fileDF <- lapply(
        tmpFiles, 
        function(x){
          data.frame(file_path = x, 
                     target = gsub(".TIF","", strsplit(basename(x), 
                                                  "S[0-9]_C[0-9]+_F[0-9]+_")[[1]][2], 
                                 fixed = T))
        }
      )
      fileDF <- do.call(rbind, fileDF)
      colnames(fileDF)[1] <- 'imgPath'
      morphFileDF <- fileDF 
      
      # per cell stats
      tmpFiles <- dir(fs::path(ProteinDir, names(fovDirs)[idx], "PerCellStats"), 
                      pattern = "S[0-9]_C[0-9]+_F[0-9]+_.*_perCell_1ChStats.csv$", recursive = F, full.names = T)
      if(length(tmpFiles)-5 != nrow(morphFileDF)){
        message(sprintf("%s: Found %d ProteinImages but %d perCell_1ChStats.", 
                        names(fovDirs)[idx], nrow(morphFileDF), length(tmpFiles)))
      }
      fileDF <- lapply(
        tmpFiles, 
        function(x){
          data.frame(file_path = x, 
                     target = gsub("_perCell_1ChStats.csv","", strsplit(basename(x), 
                                                       "S[0-9]_C[0-9]+_F[0-9]+_")[[1]][2], 
                                   fixed = T))
        }
      )
      fileDF <- do.call(rbind, fileDF)
      colnames(fileDF)[1] <- 'proteinStatsPath'
      
      # keep all images, additional stats for C902 are ignores
      morphFileDF <- merge(morphFileDF, fileDF, by = 'target', all.X = TRUE)
      
      morphFileDF[['slide']] <- slideID
      morphFileDF[['fov']] <- unname(fovDirs)[idx]
      
      return(morphFileDF)
    }
  )
  allFileDF <- do.call(rbind, allFileDF)
  
  tgrtDF <- table(allFileDF[c('fov', 'target')])
  if(sum(tgrtDF!=1)>0){
    tmpIdxDF <- which(tgrtDF !=1, arr.ind = T)
    message("All query cells must have same targets in protein image to extract them together.")
    message(sprintf("Some FOVs have missing targets: FOV %s.\n Affected targets: `%s`", 
                    paste0(rownames(tgrtDF)[unique(tmpIdxDF[, 1])], collapse = ", FOV "), 
                    paste0(colnames(tgrtDF)[unique(tmpIdxDF[, 2])], collapse = "`, `")))
  }
  
  return(allFileDF)
  
}


#' @title extractNeighborhoodImgs_fromProteinDir
#' @description extract the neighborhood images around query cells for ProteinImages
#' @param pixel_size_um pixel size in um unit of raw images
#' @param cropSize_um dimension of extracted neighborhood for query cells in um unit
#' @param cells_to_use vector of cell_ID for query cells in `c_[slide]_[fov]_[CellId]` format 
#' @param targets_to_load a vector of protein targets to load
#' @param proteinFileDF data.frame of file path for the protein images of each fov
#' @param spatColns a vector of column names for xy spatial coordinates of each cell in PerCellStats files, in pixel unit 
#' @param equalizeFOV_flag flag to equalize intensity of morph images for each channel at fov level (default = FALSE)
#' @param coreNum number of cores for mclapply parallelization 
#' @return a nested list of 2 elements, `morph`and `morphChannelIDs`, the former of which contains cell-wise image array with multi-channel for different protein targets, could be appended to the `morph` outputs of \code{extractNeighborhoodImgs} function.
#' @export 
extractNeighborhoodImgs_fromProteinDir  <- function(pixel_size_um = 0.120280945, 
                                                    cropSize_um = 50,
                                                    cells_to_use, 
                                                    targets_to_load, 
                                                    proteinFileDF, 
                                                    spatColns = c('CenterX', 'CenterY'), 
                                                    equalizeFOV_flag = FALSE, 
                                                    coreNum = 5){
  
  colns <- c('slide', 'fov', 'imgPath', 'proteinStatsPath', 'target')
  if(!all(colns %in% colnames(proteinFileDF))){
    stop(sprintf("Missing essential columns in `proteinFileDF`: `%s`.", 
                 paste0(setdiff(colns, colnames(proteinFileDF)), collapse = "`, `")))
  }
  
  if(length(spatColns) !=2){
    stop("Must provide 2 spatail column names in `spatColns`.")
  }
  
  cropSize_px <- round(cropSize_um/pixel_size_um)
  message(sprintf("Extract %d px X %d px box around each query cells.", 
                  cropSize_px+1, cropSize_px+1))
  
  
  queryDF <- lapply(
    cells_to_use, 
    function(x){
      x2 <- as.numeric(strsplit(x, '_')[[1]])
      tt <- data.frame(cell_ID= x, 
                       slide = x2[2], 
                       fov = x2[3], 
                       CellId = x2[4])
      return(tt)
    })
  queryDF <- do.call(rbind, queryDF)
  queryDF[['UID']] <- paste0(queryDF[['slide']], '_', queryDF[['fov']])
  
  # check fovs 
  proteinFileDF <- proteinFileDF[, colns, drop = F]
  proteinFileDF[['UID']] <- paste0(proteinFileDF[['slide']], '_', proteinFileDF[['fov']])
  
  sharedUIDs <- unique(intersect( queryDF[['UID']], proteinFileDF[['UID']]))
  
  if(length(sharedUIDs)<1){
    stop("Found no matched record between `cells_to_use` and `proteinFileDF`, please check if cell id is in `c_[slide]_[fov]_[CellId]` format.")
  }
  
  proteinFileDF <-  proteinFileDF[ proteinFileDF[['UID']] %in% sharedUIDs, ]
  
  # check targets shared by all query cells 
  tgrtDF <- table(proteinFileDF[c('UID', 'target')])
  if(sum(tgrtDF!=1)>0){
    tmpIdxDF <- which(tgrtDF !=1, arr.ind = T)
    message("All query cells must have same targets in protein image to extract them together.")
    message(sprintf("Some FOVs have missing targets.\nAffected UIDs: %s.\nAffected targets: `%s`", 
                    paste0(rownames(tgrtDF)[unique(tmpIdxDF[, 1])], collapse = ", "), 
                    paste0(colnames(tgrtDF)[unique(tmpIdxDF[, 2])], collapse = "`, `")))

    sharedTgrts <- colnames(tgrtDF)[setdiff(seq_len(ncol(tgrtDF)), unique(tmpIdxDF[, 2]))]
  }else{
    sharedTgrts <- colnames(tgrtDF)
  }
  missingTgrts <- setdiff(targets_to_load, sharedTgrts)
  
  if(length(missingTgrts)>0){
    stop(sprintf("%d `targets_to_load` are not shared by all query cells in `proteinFileDF`: `%s`",
                 length(missingTgrts), paste0(missingTgrts, collapse = "`, `")))
  }else{
    message(sprintf("Load protein images into multi-channel array in following order: `%s`",
                    paste0(targets_to_load, collapse = "`, `")))
  }
  
  proteinFileDF <- proteinFileDF[proteinFileDF[['target']] %in% targets_to_load, ]
  
  # get px coordinates for each cell
  tmpDF <- proteinFileDF[proteinFileDF[['target']] == targets_to_load[1], c('UID', 'proteinStatsPath'), drop = F]
  queryDF <- merge(queryDF, unique(tmpDF))

  tmpDF <- split(queryDF, as.factor(queryDF[['UID']]))
  tmpDF <- lapply(
    tmpDF, 
    function(df){
      cellStatsDF <- read.csv(df[1, 'proteinStatsPath'])
      
      if(!all(c(spatColns, 'CellId') %in% colnames(cellStatsDF))){
        stop(sprintf("Missing essential columns (`%s`) in cell stats file: `%s`.", 
                     paste0(setdiff(c(spatColns, 'CellId'), colnames(cellStatsDF)), collapse = "`, `"),
                     df[1, 'proteinStatsPath']))
      }
      
      cellStatsDF <- cellStatsDF[, c('CellId', spatColns)]
      cellStatsDF <- merge(df, cellStatsDF, by = 'CellId')
      
      return(cellStatsDF)
    }
  )
  tmpDF <- do.call(rbind, tmpDF)
  if(nrow(tmpDF) != nrow(queryDF)){
    message(spritf("%d out of %d query cells have NO cell stats information.", nrow(queryDF) - nrow(tmpDF), nrow(queryDF)))
  }
  
  queryDF <- tmpDF
  
  
  # inputs: df[, c('UID', 'CellId','cell_ID','CenterX','CenterY')], 
  # global values: cropSize_px, spatColns, equalizeFOV_flag, proteinFileDF[, c('UID', 'target','imgPath')]
  # fill 0 for crop window smaller than others 
  extractProteinPerFovFun <- function(df){
    uid <- df$UID[1]
    
    message('UID (slide_FOV) ', uid)
    
    # get protein images in same order as targets_to_load
    imgFileDF <- proteinFileDF[proteinFileDF[['UID']] == uid, ]
    rownames(imgFileDF) <- imgFileDF[['target']]
    imgFileDF <- imgFileDF[targets_to_load, ]
    
    img1 <- EBImage::readImage(imgFileDF[['imgPath']])
    fov_size <- dim(img1)[1:2]
    
    
    # get topLeft corner of crop box, clip at 1
    tmpCoord <- df[, spatColns] - floor(cropSize_px/2)
    colnames(tmpCoord) <- paste0('min_', c('x_px', 'y_px'))
    tmpCoord <- sweep(as.matrix(tmpCoord), 2, c(1,1), pmax)
    
    spatCoord <- cbind(df[, c('cell_ID', spatColns)], as.data.frame(tmpCoord))
    
    # get bottomRight corner, clip at maximum fov size 
    tmpCoord <- df[, spatColns] + ceiling(cropSize_px/2)
    colnames(tmpCoord) <- paste0('max_', c('x_px', 'y_px'))
    tmpCoord <- sweep(as.matrix(tmpCoord), 2, fov_size, pmin)
    
    spatCoord <- cbind(spatCoord, as.data.frame(tmpCoord))
    rownames(spatCoord) <- spatCoord[['cell_ID']]
    
    # equalize across channels for morph image 
    if(equalizeFOV_flag){
      img1 <- EBImage::equalize(img1, levels = 2**16)
    }
    
    # crop images for each cell
    img2 <- lapply(
      seq_len(nrow(spatCoord)), 
      function(idx){
        # return the image array by itself
        img1_crop <- img1[seq(spatCoord[idx, 'min_x_px'], 
                              spatCoord[idx, 'max_x_px']), 
                          seq(spatCoord[idx, 'min_y_px'], 
                              spatCoord[idx, 'max_y_px']), ]
        img1_crop <- img1_crop@.Data
        
        # fill 0 for region smaller than crop window 
        missX <- cropSize_px +1 - dim(img1_crop)[1]
        missY <- cropSize_px +1 - dim(img1_crop)[2]
        
        if(missX >0){
          mockData <- array(0, dim = c(missX, dim(img1_crop)[2:3]))
          img1_crop <- abind::abind(img1_crop, mockData, along = 1)
        }
        
        if(missY >0){
          mockData <- array(0, dim = c(dim(img1_crop)[1], missY, dim(img1_crop)[3]))
          img1_crop <- abind::abind(img1_crop, mockData, along = 2)
        }
        
        return(img1_crop)
      })
    names(img2) <- rownames(spatCoord)
    
    return(img2)
  }
  
  # extract info for each fov 
  perFOV_DF <- split(queryDF, as.factor(queryDF[['UID']]))
  
  
  imgData <- parallel::mclapply(perFOV_DF, 
                                mc.cores = coreNum, 
                                FUN = extractProteinPerFovFun)
  
  # get cell_id
  imgListCellIds <- unlist(lapply(imgData, function(x) names(x)))
  
  imgData <- do.call(c, imgData)
  names(imgData) <- imgListCellIds
  
  finalRes <- list(morphChannelIDs = targets_to_load, 
                   morph = imgData)
  
  return(finalRes)
} 

#' @title combineMorphList
#' @description combine 2 sets of multi-channel image array together and 
#' @param imgList1 named list of image arrays with one query cell per element, names for cell_ID
#' @param imgList2 named list of image arrays with one query cell per element, names for cell_ID
#' @return a named list of multi-channel image arrays with one query cell per element, names for cell_ID
#' @export 
combineMorphList <- function(imgList1, 
                             imgList2){
  if(is.null(names(imgList1)) || is.null(names(imgList2))){
    stop("The provided image lists must have names in cell_ID.")
  }
  sharedCells <- intersect(names(imgList1), names(imgList1))
  if(length(sharedCells)<1){
    stop("No shared cells found between two provided image lists.")
  }
  
  combinedList <- lapply(
    sharedCells, 
    function(cellName){
      img1 <- imgList1[[cellName]]
      img2 <- imgList2[[cellName]]
      
      if(!identical(dim(img1)[1:2], dim(img2)[1:2])){
        stop(sprintf("Different image dimensions between the two sets for cell: %s", cellName))
      }
      
      # if both are 1 channel 
      if(max(length(dim(img1)),  length(dim(img2)))<3){
        tmpArray <- array(0, dim = c(dim(img2)[1:2], 1))
        tmpArray[, ,1] <- img2
        img2 <- tmpArray
        rm(tmpArray)
      }
      
      # combine channel together 
      combImg <- EBImage::abind(img1, img2, along = 3)
      
      return(combImg)
    }
  )
  names(combinedList) <- sharedCells
  
  return(combinedList)
}

#' @title decodeStitchedLabelVal
#' @description decode the stitched label value in napari-CosMx stitched data to fov and CellId
#' @param labVals vector of encoded label value used in napari-CosMx stitched data 
#' @return data.frame with fov and CellId for each encoded label value
#' @export 
decodeStitchedLabelVal <- function(labVals){
  # # decode cell id and fov from label value
  # https://github.com/Nanostring-Biostats/napari-CosMx/blob/955e221c07975068a840ba6fbc4e9f68c75e7f06/src/napari_cosmx/pairing.py#L18
  labVals <- as.numeric(labVals)
  
  tmpMeta <- cbind(labVals, floor(sqrt(labVals)))
  tmpMeta <- cbind(tmpMeta, labVals - tmpMeta[, 2]**2)
  tmpMeta <- cbind(tmpMeta, tmpMeta[, 3] - tmpMeta[, 2])
  colnames(tmpMeta) <- c('labVal', 'zflr', 'fov1', 'cell2')
  
  cellMeta <- as.data.frame(tmpMeta[, c('labVal', 'zflr', 'cell2'), drop = FALSE])
  colnames(cellMeta) <- c('labVal', 'fov', 'CellId')
  
  tmpFlag <- (tmpMeta[, 'cell2'] <0)
  cellMeta[['fov']][tmpFlag] <- tmpMeta[tmpFlag, 'fov1']
  cellMeta[['CellId']][tmpFlag] <- tmpMeta[tmpFlag, 'zflr']
  
  return(cellMeta)
  
}

#' @title encodeCellId
#' @description encode fov and CellId to the stitched label value used in napari-CosMx stitched data
#' @param cellMeta data.frame of cell metadata 
#' @param cellLabel_coln column name of cell label values at FOV level in `cellMeta` (default = 'CellId')
#' @param fov_coln column name of fov in `cellMeta`
#' @return a vector of encoded label value used in napari-CosMx stitched data 
#' @export 
encodeCellId <- function(cellMeta, 
                         cellLabel_coln = 'CellId',
                         fov_coln = 'fov'){
  # # encode CellId and fov to label value in stitched data 
  # https://github.com/Nanostring-Biostats/napari-CosMx/blob/955e221c07975068a840ba6fbc4e9f68c75e7f06/src/napari_cosmx/pairing.py#L4
  
  colns <- c(cellLabel_coln, fov_coln)
  if(!all(colns %in% colnames(cellMeta))){
    stop(sprintf("Missing essential columns in `cellMeta`: `%s`.", 
                 paste0(setdiff(colns, colnames(cellMeta)), collapse = "`, `")))
  }
  
  
  tmpFlag <- (cellMeta[[fov_coln]] < cellMeta[[cellLabel_coln]])
  labVal <- cellMeta[[fov_coln]]**2 + cellMeta[[fov_coln]] + cellMeta[[cellLabel_coln]]
  labVal[tmpFlag] <- cellMeta[tmpFlag, cellLabel_coln]**2 + cellMeta[tmpFlag, fov_coln]
  
  return(labVal)
  
}


#' @title prepStitchedData  
#' @description collect information for the stitched zarr images generated by napari-CosMx plugin, link the values in labels image to the input FOV and cell_label; slower than \code{getGlobalCoords_forStitched}
#' @param napari_imgDir file path to stitched images generated by napari 
#' @param slideID slide ID that would be assigned to cells in the query flow cell
#' @return a list of 3 data.frame and 1 list
#' \itemize{
#'  \item{zarrDF}{data.frame for file path of different image channels under `napari_imgDir`}
#'  \item{fovOffsets}{data.frame for fov offset used to arrange and stitch FOVs in `napari_imgDir`}
#'  \item{cellMeta}{data.frame for cell metadata, linking original fov and CellId to the values in stitched labels image of `napari_imgDir`}
#'  \item{dimSetup}{list of 5 elements, `fov_size`,  `pixel_size_mm`, `isDASH`, `global_size`, `topLeft_mm`}
#' }
#' @export 
prepStitchedData <- function(napari_imgDir, 
                             slideID = 1){
  message("This function may take large memory and run slow due to extraction of cell_id and global coordinates directly from stitched data.\nTry `getGlobalCoords_forStitched` function if per FOV CellId and local coordinates are known.")
  napari_imgDir <- normalizePath(napari_imgDir)
  tmpFile <- fs::path(napari_imgDir, ".zattrs")
  if(!file.exists(tmpFile)){
    stop(sprintf("Cannot find `.zattrs` metadata file under `napari_imgDir`: %s", napari_imgDir))
  }
  
  studyConfigs <- jsonlite::read_json(tmpFile)[['CosMx']]

  
  # napari-CosMx flip x and y coordinates in input image 
  dimSetup <- list(fov_size = c(studyConfigs[['fov_width']], studyConfigs[['fov_height']]),
                   pixel_size_mm = ifelse(is.null(studyConfigs[['scale_mm']]), 
                                          studyConfigs[['scale_um']]/1000, 
                                          studyConfigs[['scale_mm']]))
  if(length(dimSetup[['pixel_size_mm']])!=1){
    stop("Not able to find pixel size in `.zattrs` metadata file.")
  }
  
  if(diff(dimSetup[['fov_size']])==0){
    isDASH = FALSE
  }else{
    isDASH = TRUE
  }
  
  fovOffsets <- sapply(studyConfigs[['fov_offsets']], unlist)
  rownames(fovOffsets) <- NULL  
  fovOffsets <- as.data.frame(fovOffsets[order(fovOffsets[, 'FOV']), ])  
  

  zarrDF <- data.frame(channelID = list.dirs(napari_imgDir, full.names = F, recursive = F))
  zarrDF[['file_path']] <- fs::path(napari_imgDir, zarrDF[['channelID']], "0")
  zarrDF[['type']] <- 'morph'
  zarrDF[['type']][zarrDF[['channelID']] == 'labels'] <- 'labels'
  rownames(zarrDF) <- zarrDF[['channelID']]
  
  if(any(!file.exists(zarrDF[['file_path']]))){
    stop(sprintf("Non-zoomed `0` folder is missing for channel: `%s`.", 
                 paste0(zarrDF[['channelID']][!file.exists(zarrDF[['file_path']])], 
                        collapse = "`, `")))
  }  
  
  
  # base location 
  if(isDASH){
    topLeft_mm <- c(-max(fovOffsets[['X_mm']]), min(fovOffsets[['Y_mm']]))
  }else{
    topLeft_mm <- c(min(fovOffsets[['Y_mm']]), -max(fovOffsets[['X_mm']]))
  }
  
  dimSetup[['isDASH']] <- isDASH
  dimSetup[['topLeft_mm']] <- topLeft_mm
  
  # get cell id and coordinates
  labArray <- Rarr::zarr_overview(zarrDF['labels', 'file_path'], as_data_frame = T)
  # y, x
  globalSize <- labArray$dim[[1]]
  dimSetup[['global_size']] <- globalSize 
  
  # read in label at bin =2, need high memory 
  labArray <- Rarr::read_zarr_array(fs::path(dirname(zarrDF['labels', 'file_path']), '1'), 
                                    index = list(seq_len(floor(globalSize[1]/2)), 
                                                 seq_len(floor(globalSize[2]/2))))
  allCells <- setdiff(unique(labArray), 0)

  # # decode cell id and fov from label value
  # https://github.com/Nanostring-Biostats/napari-CosMx/blob/955e221c07975068a840ba6fbc4e9f68c75e7f06/src/napari_cosmx/pairing.py#L18
  cellMeta <- decodeStitchedLabelVal(allCells)
  
  cellMeta[['slide']] <- slideID
  cellMeta[['cell_ID']] <- paste0('c_', slideID, '_', cellMeta[['fov']], '_', cellMeta[['CellId']])
  
  # get centroid coordinates for each cell
  tmpMeta <- sapply(
    allCells, 
    function(cl){
      df <- which(labArray == cl, arr.ind = T)
      df <- apply(df, 2, median)
      return(df)
    }
  )
  # get coords for bin =1 case
  cellMeta <- cbind(tmpMeta*2, cellMeta)
  colnames(cellMeta)[1:2] <- c('global_y_px', 'global_x_px')
  
  res <- list(zarrDF = zarrDF, 
              fovOffsets = fovOffsets, 
              cellMeta = cellMeta, 
              dimSetup = dimSetup)
  return(res)
}


#' @title getGlobalCoords_forStitched  
#' @description convert local coordinates at FOV level into global coordinates and link cell_ID to label values used in napari-CosMx stitched data
#' @param napari_imgDir file path to stitched images generated by napari 
#' @param cellMeta data.frame for cell metadata of all cells in napari-CosMx stitched data, can be combined info of all Cell_Stats files in CellStatsDir.
#' @param cellID_coln column name of cell_ID in `cellMeta`
#' @param cellLabel_coln column name of cell label values at FOV level in `cellMeta` (default = 'CellId')
#' @param fov_coln column name of fov in `cellMeta`
#' @param spatColns a vector of column names for local xy spaitial coordinates of each cell in `cellMeta`, in pixel unit 
#' @return a list of 3 data.frame and 1 list
#' \itemize{
#'  \item{zarrDF}{data.frame for file path of different image channels under `napari_imgDir`}
#'  \item{fovOffsets}{data.frame for fov offset used to arrange and stitch FOVs in `napari_imgDir`}
#'  \item{cellMeta}{data.frame for cell metadata, linking original fov and CellId to the values in stitched labels image of `napari_imgDir`}
#'  \item{dimSetup}{list of 5 elements, `fov_size`,  `pixel_size_mm`, `isDASH`, `global_size`, `topLeft_mm`}
#' }
#' @export 
getGlobalCoords_forStitched <- function(napari_imgDir, 
                                        cellMeta, 
                                        cellID_coln = 'cell_ID', 
                                        cellLabel_coln = 'CellId',
                                        fov_coln = 'fov',
                                        spatColns = c('CenterX', 'CenterY')){
  colns <- c(cellID_coln, cellLabel_coln, fov_coln, spatColns)
  if(!all(colns %in% colnames(cellMeta))){
    stop(sprintf("Missing essential columns in `cellMeta`: `%s`.", 
                 paste0(setdiff(colns, colnames(cellMeta)), collapse = "`, `")))
  }
  
  if(length(spatColns) !=2){
    stop("Must provide 2 spatail column names in `spatColns`.")
  }
  # drop unnecessary columns
  cellMeta <- cellMeta[, colns, drop = F]
  
  # encode CellId and fov to label value in stitched data 
  cellMeta[['labVal']] <- encodeCellId(cellMeta = cellMeta, 
                                       cellLabel_coln = cellLabel_coln, 
                                       fov_coln = fov_coln)
  
  # check stitched data folder and get dimension setup
  napari_imgDir <- normalizePath(napari_imgDir)
  tmpFile <- fs::path(napari_imgDir, ".zattrs")
  if(!file.exists(tmpFile)){
    stop(sprintf("Cannot find `.zattrs` metadata file under `napari_imgDir`: %s", napari_imgDir))
  }
  
  studyConfigs <- jsonlite::read_json(tmpFile)[['CosMx']]
  
  
  # napari-CosMx flip x and y coordinates in input image 
  dimSetup <- list(fov_size = c(studyConfigs[['fov_width']], studyConfigs[['fov_height']]),
                   pixel_size_mm = ifelse(is.null(studyConfigs[['scale_mm']]), 
                                          studyConfigs[['scale_um']]/1000, 
                                          studyConfigs[['scale_mm']]))
  if(length(dimSetup[['pixel_size_mm']])!=1){
    stop("Not able to find pixel size in `.zattrs` metadata file.")
  }
  
  
  fovOffsets <- sapply(studyConfigs[['fov_offsets']], unlist)
  rownames(fovOffsets) <- NULL  
  fovOffsets <- as.data.frame(fovOffsets[order(fovOffsets[, 'FOV']), ])  
  
  
  zarrDF <- data.frame(channelID = list.dirs(napari_imgDir, full.names = F, recursive = F))
  zarrDF[['file_path']] <- fs::path(napari_imgDir, zarrDF[['channelID']], "0")
  zarrDF[['type']] <- 'morph'
  zarrDF[['type']][zarrDF[['channelID']]== 'labels'] <- 'labels'
  rownames(zarrDF) <- zarrDF[['channelID']]
  
  if(any(!file.exists(zarrDF[['file_path']]))){
    stop(sprintf("Non-zoomed `0` folder is missing for channel: `%s`.", 
                 paste0(zarrDF[['channelID']][!file.exists(zarrDF[['file_path']])], 
                        collapse = "`, `")))
  }  
  
  labArray <- Rarr::zarr_overview(zarrDF['labels', 'file_path'], as_data_frame = T)
  # y, x
  dimSetup[['global_size']] <- labArray$dim[[1]]
  
  
  # get global coordinates 
  if(diff(dimSetup[['fov_size']])==0){
    dimSetup[['isDASH']] <- FALSE
    topLeft_mm <- c(min(fovOffsets[['Y_mm']]), -max(fovOffsets[['X_mm']]))
    spatCoords <- fovOffsets[, c('FOV', 'X_mm', 'Y_mm')]
    
  }else{
    dimSetup[['isDASH']] <- TRUE
    topLeft_mm <- c(-max(fovOffsets[['X_mm']]), min(fovOffsets[['Y_mm']]))
    spatCoords <- fovOffsets[, c('FOV', 'X_mm', 'Y_mm')]
    # flip stage axis and orientation for DASH
    spatCoords[c(1,2)] <- (- spatCoords[c(2,1)] )
  }
  colnames(spatCoords) <- c(fov_coln, 'x_stage_mm', 'y_stage_mm')
  dimSetup[['topLeft_mm']] <- topLeft_mm
  
  spatCoords <- merge(cellMeta[, c(cellID_coln, fov_coln, spatColns)], 
                      spatCoords, by = fov_coln)
  spatCoords[['global_x_px']] <- spatCoords[[spatColns[1]]] - spatCoords[['x_stage_mm']]/dimSetup[['pixel_size_mm']] - topLeft_mm[2]/dimSetup[['pixel_size_mm']] 
  spatCoords[['global_y_px']] <- spatCoords[[spatColns[2]]] + spatCoords[['y_stage_mm']]/dimSetup[['pixel_size_mm']] - topLeft_mm[1]/dimSetup[['pixel_size_mm']] 
  
  spatCoords[['global_x_mm']] <- spatCoords[['global_x_px']]*dimSetup[['pixel_size_mm']] 
  spatCoords[['global_y_mm']] <- spatCoords[['global_y_px']]*dimSetup[['pixel_size_mm']] 
  rownames(spatCoords) <- spatCoords[[cellID_coln]]
  
  cellMeta <- cbind(cellMeta, 
                    spatCoords[cellMeta[[cellID_coln]], 
                               c('global_x_px', 'global_y_px', 'global_x_mm', 'global_y_mm')])
  
  res <- list(zarrDF = zarrDF, 
              fovOffsets = fovOffsets, 
              cellMeta = cellMeta, 
              dimSetup = dimSetup)
  return(res)
}


#' @title extractNeighborhoodImgs_fromStitched
#' @description extract the neighborhood images around query cells from napari-CosMx stitched data 
#' @param pixel_size_um pixel size in um unit of raw images
#' @param cropSize_um dimension of extracted neighborhood for query cells in um unit
#' @param cells_to_use vector of cell_ID for query cells in `c_[slide]_[fov]_[CellId]` format 
#' @param zarrDF data.frame of file path for the paired napari-CosMx stitched data
#' @param global_size a vector for the dimension of stitched image in pixel unit and c(column/width, rows/height) order, output of \code{prepStitchedData} or \code{getGlobalCoords_forStitched} function
#' @param cellMeta data.frame for cell metadata of cells in napari-CosMx stitched data, output of \code{prepStitchedData} or \code{getGlobalCoords_forStitched} function
#' @param cellID_coln column name of cell_ID in `cellMeta`
#' @param cellLabel_coln column name of cell label values in `cellMeta` (default = 'labVal')
#' @param spatColns a vector of column names for global xy spatial coordinates of each cell in `cellMeta`, in pixel unit and c(column index, row index) order
#' @param flip_axes flag to flip between row and column axes of extracted image upon return (default = TRUE)
#' @param equalize_flag flag to equalize intensity of morph images for each channel at each query cell level (default = FALSE)
#' @param coreNum number of cores for mclapply parallelization 
#' @return a nested list of 3 elements, `labels`, `morph`, `morphChannelIDs`, the former of which contain list of cell-wise image array, the 3rd element is a vector of channelID in `morph`
#' @export 
extractNeighborhoodImgs_fromStitched <- function(pixel_size_um = 0.120280945, 
                                                 cropSize_um = 50,
                                                 cells_to_use, 
                                                 zarrDF, 
                                                 global_size, 
                                                 cellMeta, 
                                                 cellID_coln = 'cell_ID', 
                                                 cellLabel_coln = 'labVal',
                                                 spatColns = c('global_x_px', 'global_y_px'), 
                                                 flip_axes = TRUE, 
                                                 equalize_flag = FALSE,
                                                 coreNum = 5){
  colns <- c(cellID_coln, cellLabel_coln, spatColns)
  if(!all(colns %in% colnames(cellMeta))){
    stop(sprintf("Missing essential columns in `cellMeta`: `%s`.", 
                 paste0(setdiff(colns, colnames(cellMeta)), collapse = "`, `")))
  }
  cellMeta <- cellMeta[cellMeta[[cellID_coln]] %in% cells_to_use, colns, drop = F]
  
  if(nrow(cellMeta)<1){
    stop(sprintf("Cannot find `cells_to_use` in the cell_ID column `%s` of provided `cellMeta`.", cellID_coln))
  }else{
    message(sprintf("Found %d `cells_to_use` in the provided `cellMeta`.", nrow(cellMeta)))
  }
  
  if(length(spatColns) !=2){
    stop("Must provide 2 spatail column names in `spatColns`.")
  }
  
  colns <- c('channelID', 'file_path', 'type')
  if(!all(colns %in% colnames(zarrDF))){
    stop(sprintf("Missing essential columns in `zarrDF`: `%s`.", 
                 paste0(setdiff(colns, colnames(zarrDF)), collapse = "`, `")))
  }
  
  zarrDF <- zarrDF[file.exists(zarrDF[['file_path']]), ]
  if(nrow(zarrDF)<1){
    stop("File not exists for the provided `zarrDF[['file_path']]`. ")
  }
  
  if(sum(zarrDF[['type']] == 'labels')!=1){
    stop(sprintf("Found %d labels channel in the provided `zarrDF[['file_path']]`. ",
                 sum(zarrDF[['type']] == 'labels')))
  }

  zarrDF[['chanOrder']] <- 0
  zarrDF[['chanOrder']][zarrDF[['type']] == 'morph'] <- seq_len(sum(zarrDF[['type']] == 'morph'))
  zarrDF <- zarrDF[order(zarrDF[['chanOrder']]), ]
  message(sprintf("Extract images for channels: `%s`.", 
                  paste0(zarrDF[['channelID']], collapse = "`, `")))
  rownames(zarrDF) <- zarrDF[['channelID']]
  
  cropSize_px <- round(cropSize_um/pixel_size_um)
  message(sprintf("Extract %d px X %d px box around each query cells.", 
                  cropSize_px+1, cropSize_px+1))
  
  
  # get topLeft corner of crop box, clip at 1
  tmpCoord <- round(cellMeta[, spatColns]) - floor(cropSize_px/2)
  colnames(tmpCoord) <- paste0('min_', c('x_px', 'y_px'))
  tmpCoord <- sweep(as.matrix(tmpCoord), 2, c(1,1), pmax)
  
  cellMeta <- cbind(cellMeta, as.data.frame(tmpCoord))
  
  # get bottomRight corner, clip at maximum fov size 
  tmpCoord <- round(cellMeta[, spatColns]) + ceiling(cropSize_px/2)
  colnames(tmpCoord) <- paste0('max_', c('x_px', 'y_px'))
  tmpCoord <- sweep(as.matrix(tmpCoord), 2, global_size[2:1], pmin)
  
  cellMeta <- cbind(cellMeta, as.data.frame(tmpCoord))
  rownames(cellMeta) <- cellMeta[[cellID_coln]]
  
  # inputs: cellMeta[, c(cellID_coln, cellLabel_coln, spatColns,'min_x_px','min_y_px', 'max_x_px', 'min_y_px')], 
  # global values: cropSize_px, spatColns, equalize_flag, flip_axes, cellID_coln, cellLabel_coln, zarrDF
  # fill 0 for crop window smaller than others 
  extractPerCellFun <- function(cMeta){
    message('cell ', cMeta[[cellID_coln]])
    
    # get images 
    img1_crop <- lapply(
      zarrDF[['channelID']], 
      function(ch){
        cellZarr <- Rarr::read_zarr_array(
          zarrDF[ch, 'file_path'], 
          index = list(seq(cMeta[['min_y_px']],
                           cMeta[['max_y_px']]),
                       seq(cMeta[['min_x_px']],
                           cMeta[['max_x_px']])))
        if(zarrDF[ch, 'type'] != 'labels'){
          # not label, normalized the intensity to 16-bit
          cellZarr <- cellZarr/2**16
        }
        return(cellZarr)
      }
    )
    
    if(length(img1_crop)>1){
      img1_crop <- abind::abind(img1_crop, along = 0)
      img1_crop <- aperm(img1_crop, c(2,3,1))
    }else{
      img1_crop <- img1_crop[[1]]
      tmpArray <- array(0, dim = c(dim(img1_crop), 1))
      tmpArray[, , 1] <- img1_crop
      img1_crop <- tmpArray
      rm(tmpArray)
    }
    
    
    # fill 0 for region smaller than crop window 
    missX <- cropSize_px +1 - dim(img1_crop)[1]
    missY <- cropSize_px +1 - dim(img1_crop)[2]
    
    if(missX >0){
      mockData <- array(0, dim = c(missX, dim(img1_crop)[2:3]))
      img1_crop <- abind::abind(img1_crop, mockData, along = 1)
    }
    
    if(missY >0){
      mockData <- array(0, dim = c(dim(img1_crop)[1], missY, dim(img1_crop)[3]))
      img1_crop <- abind::abind(img1_crop, mockData, along = 2)
    }
    
    # flip xy axes of image array
    if(flip_axes){
      img1_crop <- aperm(img1_crop, c(2,1,3))
    }
      
    # process morph
    if(dim(img1_crop)[3]>1){
      imgMorph <- img1_crop[, , seq(2, dim(img1_crop)[3])]
      
      # equalize across channels for morph image 
      if(equalize_flag){
        imgMorph <- EBImage::equalize(imgMorph, levels = 2**16)
      }
 
    }else{
      imgMorph <- NA
    }
    
    cellRes <- list(labels = img1_crop[, , 1], 
                    morph = imgMorph)

    return(cellRes)
  }
  
  # extract info for each cell 
  tmpDF <- split(cellMeta, as.factor(cellMeta[[cellID_coln]]))
  imgData <- parallel::mclapply(tmpDF, 
                                mc.cores = coreNum, 
                                FUN = extractPerCellFun)
  
  names(imgData) <- names(tmpDF)

  # restructure to return 2 list for each element separately 
  finalRes <- list(labels = lapply(imgData, '[[', 'labels'), 
                   morph = lapply(imgData, '[[', 'morph'),
                   morphChannelIDs = zarrDF[['channelID']][-1])
  
  return(finalRes)
}


  