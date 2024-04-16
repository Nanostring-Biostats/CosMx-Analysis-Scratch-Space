#' @title createComposite_singleFrame
#' @description create composite RGB images from the input single-frame multi-channel image array 
#' @param img array or EBImage class for single-frame image
#' @param vizChannels vector of index of channels to show
#' @param channelColors vector of monochrome to use for each channel in vizChannels
#' @param equalizeFlags vector of logic values to equalize each channel in vizChannels
#' @return a rgb image array of 3 channels in EBImage class
#' @export 
createComposite_singleFrame <- function(img,
                                        vizChannels = NULL, 
                                        channelColors = NULL, 
                                        equalizeFlags = NULL ){
  img <- EBImage::Image(img, colormode = EBImage::Grayscale)
  img.dim <- dim(img)
  
  
  if(length(img.dim)<3){
    # single channel, single frame 
    tmpImg <- array(0, dim = c(img.dim, 1))
    tmpImg[, , 1] <- img
    img <- EBImage::Image(tmpImg)
    rm(tmpImg)
  }else if(length(img.dim)>3){
    # multi-channel, multi frame 
    stop("`img` must be a signle-frame image. ")
  }
  
  img.dim <- dim(img)
  chanNum <- img.dim[3]
  
  
  if(!all(vizChannels %in% seq_len(chanNum))){
    stop(sprintf("The provided `vizChannels` has channel index not matched with %d channels in input `img`: %s. ", 
                 chanNum, paste0(vizChannels, collapse = ", ")))
  }
  
  if(length(vizChannels) != length(channelColors)){
    stop("The `channelColors` provided must have same length as `vizChannels`.")
  }
  
  if(!is.null(equalizeFlags)){
    if(length(vizChannels) != length(equalizeFlags)){
      stop("The `equalizeFlag` provided must have same length as `vizChannels`.")
    }
    
    if(!is.logical(equalizeFlags)){
      stop("The `equalizeFlag` provided must be a vector of logic flag for channels in `vizChannels`.")
    }
    
  }
  
  colorCode <- matrix(0, nrow = 3, ncol = 7, 
                      dimnames = list(c("asred", "asgreen", "asblue"), 
                                      c('Red', 'Green', 'Blue', 'Cyan', 'Magenta', 'Yellow', 'Greys')))
  colorCode['asred', c('Red', 'Magenta', 'Yellow', 'Greys')] <- 1
  colorCode['asgreen', c('Green', 'Cyan', 'Yellow', 'Greys')] <- 1
  colorCode['asblue', c('Blue', 'Cyan', 'Magenta', 'Greys')] <- 1
  
  if(!all(channelColors %in% colnames(colorCode))){
    stop(sprintf("The `channelColors` must choose from `%s`.", 
                 paste0(colnames(colorCode), collapse = "', '")))
  }
  
  # convert each channel into rgb
  rgbImgs <- lapply(
    seq_along(vizChannels), 
    function(idx){
      vizIdx <- vizChannels[idx]
      vizColorC <- colorCode[, channelColors[idx]]
      vizColorC <- which(vizColorC ==1)
      
      toEqulize <- equalizeFlags[idx]
      if(is.null(toEqulize)){
        toEqulize <- FALSE
      }
      
      tmpImg <- array(0, dim = c(img.dim[1:2], 3))
      
      if(toEqulize){
        tmpImg[, , vizColorC] <- EBImage::equalize(img[, , vizIdx], levels = 2**16)
      }else{
        tmpImg[, , vizColorC] <- img[, , vizIdx]
      }
      
      return(tmpImg)
    }
  )
  
  # flatten all vizChannels together 
  if(length(rgbImgs)>1){
    rgbImgs <- Reduce(function(x,y) {x+y}, rgbImgs)
  }else{
    rgbImgs <- rgbImgs[[1]]
  }
  
  rgbImgs <- EBImage::Image(rgbImgs, colormode = EBImage::Color)
  
  return(rgbImgs)
}


#' @title createComposite_multiFrame
#' @description create composite RGB images from the input single-frame multi-channel image array 
#' @param imgObj array or EBImage class for multi-frame image
#' @param vizChannels vector of index of channels to show
#' @param channelColors vector of monochrome to use for each channel in vizChannels
#' @param equalizeFlags vector of logic values to equalize each channel in vizChannels
#' @return a multi-frame rgb image array of 3 channels in EBImage class 
#' @export 
createComposite_multiFrame <- function(imgObj,
                                       vizChannels = NULL, 
                                       channelColors = NULL, 
                                       equalizeFlags = NULL ){
  obj.dim <- dim(imgObj)
  
  if(length(obj.dim) <3){
    stop("The provided `imgObj` has only 2 dimensions, use `createComposite_singleFrame()` function instead.")
  } else if(length(obj.dim) ==3){
    message("The provided `imgObj` has only 3 dimensions, assuming it's multi-frame images of only 1 channel.")
    tmpArray <- array(0, dim = c(obj.dim, 1))
    tmpArray[, , , 1]  <- imgObj
    imgObj <- aperm(tmpArray, c(1,2,4,3))
    rm(tmpArray)
  }else if(length(obj.dim) > 4){
    stop("The provided `imgObj` has more than 4 dimensions, beyond the accepted range.")
  }
  
  obj.dim <- dim(imgObj)
  
  comMorphImg <- lapply(
    seq_len(obj.dim[4]), 
    function(idx){
      img <- imgObj[, , , idx]
      rgbImg <- createComposite_singleFrame(
        img = img,
        vizChannels = vizChannels, 
        channelColors = channelColors, 
        equalizeFlags = equalizeFlags)
      return(rgbImg)
    }
  )
  
  if(obj.dim[4]>1){
    comMorphImg <- abind::abind(comMorphImg, along = 0)
    comMorphImg <- aperm(comMorphImg, c(2,3,4,1))
  }else{
    comMorphImg <- comMorphImg[[1]]
    tmpArray <- array(0, dim = c(dim(comMorphImg), 1))
    tmpArray[, , , 1] <- comMorphImg
    comMorphImg <- tmpArray
    rm(tmpArray)
  }
  
  comMorphImg <- EBImage::Image(comMorphImg, colormode = EBImage::Color)
  
  return(comMorphImg)
}


#' @title combineImgByGroup
#' @description convert list of image arrays into multi-frame image object, one object per group 
#' @param imgList named list of image arrays with one query cell per element, names for cell_ID
#' @param groupVector named vector of group ID with names to be cell_ID, values to be group ID. 
#' @param imgType type of input image in `imgList`, labels for 2D matrix, morph for multi-channel 3D array
#' @return a nested list of 2 elements, `cell_IDs` and `combImgObj`, whose elements are organized by group
#' @export 
combineImgByGroup <- function(imgList, 
                              groupVector, 
                              imgType = c('labels', 'morph')){
  imgType <- match.arg(imgType, c('labels', 'morph'))
  
  if(length(groupVector) != length(imgList)){
    stop("Must have same length in `imgList`, `groupVector`, with each element for one query cell.")
  }
  
  cells_to_use <- names(imgList)
  if(is.null(cells_to_use)){
    cells_to_use <- names(groupVector)
    if(is.null(cells_to_use)){
      stop("No cell_ID provided in the names of either `imgList` or `groupVector`.")
    }else{
      names(imgList) <- cells_to_use
    }
  }else{
    if(is.null(names(groupVector))){
      names(groupVector) <- cells_to_use
      message("No names in provided `groupVector`, use the names of `imgList` as cell_ID.")
    }else{
      if(!all(names(groupVector) %in% cells_to_use)){
        stop("The names of `imgList` and `groupVector` are not matched.")
      }
    }
  }
  
  # combined cells within same group
  cellList <- split(names(groupVector), as.factor(groupVector))
  
  
  imgData <- lapply(
    names(cellList), 
    function(cl){
      cells_to_plot <- cellList[[cl]]
      
      # combine all cell label image
      if(imgType == 'labels'){
        
        combImg <- lapply(
          imgList[cells_to_plot], 
          function(img){
            if(length(dim(img))>2){
              message("More than 1 channel provided for labels, use only the 1st channel")
              img <- img[, , 1]
            }
            colorImg <- EBImage::Image(img)
            return(colorImg)
          })
        if(length(cells_to_plot)>1){
          combImg <- Reduce(function(x,y) EBImage::combine(x,y), combImg)
        }else{
          combImg <- combImg[[1]]
          tmpArray <- array(0, dim = c(dim(combImg), 1))
          tmpArray[, , 1] <- combImg
          combImg <- EBImage::Image(tmpArray)
          rm(tmpArray)
        }
      }
      
      
      # combine all rgb morph image 
      if(imgType == 'morph'){
        combImg <- lapply(
          imgList[cells_to_plot], 
          function(img){
            if(length(dim(img))>3){
              message("More than 1 frame provided for morph, use only the 1st frame")
              img <- img[, , ,1]
            }
            colorImg <- EBImage::Image(img)
            return(colorImg)
          })
        if(length(cells_to_plot)>1){
          combImg <- abind::abind(combImg, along = 0)
          combImg <- aperm(combImg, c(2,3,4,1))
          combImg <- EBImage::Image(combImg)
        }else{
          combImg <- combImg[[1]]
          tmpArray <- array(0, dim = c(dim(combImg), 1))
          tmpArray[, , , 1] <- combImg
          combImg <- EBImage::Image(tmpArray)
          rm(tmpArray)
        }
        
      }
      
      return(combImg)
      
    }
  )
  
  names(imgData) <- names(cellList)
  
  finalRes <- list(cell_IDs = cellList, 
                   combImgObj = imgData)
  
  return(finalRes)
  
} 



#' @title overlayCellBorder_imgList
#' @description overlay cell borders on top of RGB images which could be outputs of \code{\link{createComposite_singleFrame}}, each input image has only 1 frame
#' @param rgbMorphImgs list of rgb image arrays with one query cell per element, each rgb image array has at least 3 channels with value in range of 0-1 for RGB color
#' @param labelsImgs list of cell label image arrays with one query cell per element, each label image array is a 2D matrix with value in integer indicating cell id
#' @param groupVector named vector of group ID with names to be cell_ID, values to be group ID. 
#' @param morphToShow a vector of 3 elements for the channel index of each RGB channel in `rgbMorphImgs`.
#' @param border_color hex code for color of cell borders on overlay (default = '#ffffff')
#' @param export_plots a vector of different output images that would be exported to png file on disk 
#' @param plotresultdir output folder for plot results 
#' @return a nested list with 4 elements, `cell_IDs`, `combOverlays`, `combLabels`, `combRGBs`, which has one group per element.   
#' @export 
overlayCellBorder_imgList <- function(rgbMorphImgs, 
                                      labelsImgs, 
                                      groupVector, 
                                      morphToShow = c(1,2,3),
                                      border_color =  '#ffffff', 
                                      export_plots = c('none', 'morph', 'labels', 'overlay'),
                                      plotresultdir = getwd()){
  
  export_plots <- intersect(export_plots, c('none', 'morph', 'labels', 'overlay'))
  if(length(export_plots)<1){
    stop("`explot_plots` must be either `none` or elements in c('morph', 'labels', 'overlay').")
  }
  
  if('none' %in% export_plots){
    explotFlag <- FALSE
    export_plots <- 'none'
  }else{
    explotFlag <- TRUE
  }
  
  
  cellNum <- unique(c(length(groupVector), 
                      length(rgbMorphImgs), 
                      length(labelsImgs)))
  if(length(cellNum)>1){
    stop("Must have same length in `rgbMorphImgs`, `labelsImgs`, `groupVector`, with each element for one query cell.")
  }
  
  cells_to_use <- names(rgbMorphImgs)
  if(is.null(cells_to_use)){
    cells_to_use <- names(groupVector)
    if(is.null(cells_to_use)){
      stop("No cell_ID provided in the names of either `rgbMorphImgs` or `groupVector`.")
    }else{
      names(rgbMorphImgs) <- cells_to_use
    }
  }else{
    if(is.null(names(groupVector))){
      names(groupVector) <- cells_to_use
      message("No names in provided `groupVector`, use the names of `rgbMorphImgs` as cell_ID.")
    }else{
      if(!all(names(groupVector) %in% cells_to_use)){
        stop("The names of `rgbMorphImgs` and `groupVector` are not matched.")
      }
    }
  }
  
  
  if(is.null(names(labelsImgs))){
    names(labelsImgs) <- cells_to_use
  }else{
    if(!all(names(labelsImgs) %in% cells_to_use)){
      stop("The names of `labelsImgs` and `groupVector` are not matched.")
    }
  }
  
  if(explotFlag){
    if(!dir.exists(plotresultdir)) dir.create(plotresultdir)
    message(sprintf('Export combined plots per group for `%s` images.', 
                    paste0(export_plots, collapse = "`, `")))
  }
  
  if(length(morphToShow)!=3){
    stop("`morphToShow` must has 3 elements for channel index of R, B, G channels, respectively.")
  }
  
  
  # combined cells within same group
  cellList <- split(names(groupVector), as.factor(groupVector))
  
  
  imgData <- lapply(
    names(cellList), 
    function(cl){
      cells_to_plot <- cellList[[cl]]
      
      # combine all cell label image
      comLabImg <- lapply(
        labelsImgs[cells_to_plot], 
        function(img){
          if(length(dim(img))>2){
            message("More than 1 channel provided for labels, use only the 1st channel")
            img <- img[, , 1]
          }
          colorImg <- EBImage::Image(img)
          return(colorImg)
        })
      if(length(cells_to_plot)>1){
        comLabImg <- Reduce(function(x,y) EBImage::combine(x,y), comLabImg)
      }else{
        comLabImg <- comLabImg[[1]]
        tmpArray <- array(0, dim = c(dim(comLabImg), 1))
        tmpArray[, , 1] <- comLabImg
        comLabImg <- EBImage::Image(tmpArray)
        rm(tmpArray)
      }
      
      
      
      # combine all rgb morph image 
      comMorphImg <- lapply(
        rgbMorphImgs[cells_to_plot], 
        function(img){
          if(length(dim(img))>3){
            message("More than 1 frame provided for morph, use only the 1st frame")
            img <- img[, , ,1]
          }
          colorImg <- EBImage::Image(img)
          return(colorImg)
        })
      if(length(cells_to_plot)>1){
        comMorphImg <- abind::abind(comMorphImg, along = 0)
        comMorphImg <- aperm(comMorphImg, c(2,3,4,1))
        comMorphImg <- EBImage::Image(comMorphImg)
      }else{
        comMorphImg <- comMorphImg[[1]]
        tmpArray <- array(0, dim = c(dim(comMorphImg), 1))
        tmpArray[, , , 1] <- comMorphImg
        comMorphImg <- EBImage::Image(tmpArray)
        rm(tmpArray)
      }
      
      
      EBImage::colorMode(comMorphImg) <- EBImage::Color
      
      # paint cell borders on RGB morphology 
      if(!all(morphToShow %in% seq_len(dim(comMorphImg)[3]))){
        stop("The input `morphToShow` is beyond the channel range of `rgbMorphImgs`.")
      }
      segImg <-  EBImage::paintObjects(comLabImg, comMorphImg[, , morphToShow, ], col= border_color)
      
      if(explotFlag){
        imagesPerRow <- round(sqrt(length(cells_to_plot)))
        
        # Get dimensions of combined image
        img.dim <- dim(comLabImg)
        
        xcenter <- img.dim[1]/2 + ((seq_len(img.dim[3])-1)%%imagesPerRow)*img.dim[1]
        yborder <- ceiling(seq_len(img.dim[3])/imagesPerRow-1)*img.dim[2]
        
        if('labels' %in% export_plots){
          png(fs::path(plotresultdir, paste0('group_', cl, '__cellLabel.png')), 
              width = 3 * imagesPerRow, height = 3*ceiling(length(cells_to_plot)/imagesPerRow), 
              units = "in", res = 300)
          
          EBImage::display(EBImage::colorLabels(comLabImg), all = TRUE, 
                           method = 'raster', nx = imagesPerRow)
          rect(xcenter - img.dim[1]*0.4, yborder+img.dim[2]*0.2,
               xcenter + img.dim[1]*0.4, yborder,
               col = "#FFFFFF80", border = NA)
          text(x = xcenter, y = yborder, 
               labels = cells_to_plot, adj = c(1, 0.5), pos = 1, col = "black")
          
          dev.off()
        }
        
        if('morph' %in% export_plots){
          png(fs::path(plotresultdir, paste0('group_', cl, '__morphRGB_', 
                                             paste0(morphToShow, collapse = "-"), 
                                             '.png')), 
              width = 3 * imagesPerRow, height = 3*ceiling(length(cells_to_plot)/imagesPerRow), 
              units = "in", res = 300)
          
          EBImage::display(comMorphImg[, , morphToShow, ], method = 'raster', 
                           all = T, nx = imagesPerRow)
          rect(xcenter - img.dim[1]*0.4, yborder+img.dim[2]*0.2,
               xcenter + img.dim[1]*0.4, yborder,
               col = "#FFFFFF80", border = NA)
          text(x = xcenter, y = yborder, 
               labels = cells_to_plot, adj = c(1, 0.5), pos = 1, col = "black")
          
          dev.off()
          
        }
        
        if('overlay' %in% export_plots){
          png(fs::path(plotresultdir, paste0('group_', cl, '__border_on_morph.png')), 
              width = 3 * imagesPerRow, height = 3*ceiling(length(cells_to_plot)/imagesPerRow), 
              units = "in", res = 300)
          
          EBImage::display(segImg, all = TRUE, 
                           method = 'raster', nx = imagesPerRow)
          rect(xcenter - img.dim[1]*0.4, yborder+img.dim[2]*0.2,
               xcenter + img.dim[1]*0.4, yborder,
               col = "#FFFFFF80", border = NA)
          text(x = xcenter, y = yborder, 
               labels = cells_to_plot, adj = c(1, 0.5), pos = 1, col = "black")
          
          dev.off()
        }
      }
      
      
      return(list(combOverlays = segImg,
                  combLabels = comLabImg, 
                  combRGBs = comMorphImg))
      
    }
  )
  
  names(imgData) <- names(cellList)
  
  finalRes <- list(cell_IDs = cellList, 
                   combOverlays = lapply(imgData, '[[', 'combOverlays'),
                   combLabels = lapply(imgData, '[[', 'combLabels'),
                   combRGBs = lapply(imgData, '[[', 'combRGBs'))
  
  return(finalRes)
}



#' @title overlayCellBorder_imgObj
#' @description overlay cell borders on top of RGB images which could be outputs of \code{\link{createComposite_multiFrame}}, each input imageObj has multiple frames 
#' @param rgbMorphObj multi-frame rgb image array with one query cell per frame, each frame has at least 3 channels with values in range of 0-1 for RGB color
#' @param labelsObj multi-frame cell label image array with one query cell per frame, each frame is a 2D matrix with value in integer indicating cell id
#' @param cell_ids a vector of cell_id for each frame in input image obj (default = NULL)
#' @param morphToShow a vector of 3 elements for the channel index of each RGB channel in `rgbMorphImgs`.
#' @param border_color hex code for color of cell borders on overlay (default = '#FFFFFF')
#' @param export_plots a vector of different output images that would be exported to png file on disk 
#' @param plotresultdir output folder for plot results 
#' @param fileprefix file name prefix for output pngs 
#' @return a multi-frame rgb image array with cell border overlay on top of morph  
#' @export 
overlayCellBorder_imgObj<- function(rgbMorphObj, 
                                    labelsObj, 
                                    cell_ids = NULL, 
                                    morphToShow = c(1,2,3),
                                    border_color =  '#FFFFFF', 
                                    export_plots = c('none', 'morph', 'labels', 'overlay'),
                                    plotresultdir = getwd(), 
                                    fileprefix = "group1"){
  
  export_plots <- intersect(export_plots, c('none', 'morph', 'labels', 'overlay'))
  if(length(export_plots)<1){
    stop("`explot_plots` must be either `none` or elements in c('morph', 'labels', 'overlay').")
  }
  
  if('none' %in% export_plots){
    explotFlag <- FALSE
    export_plots <- 'none'
  }else{
    explotFlag <- TRUE
  }
  
  if(explotFlag){
    if(!dir.exists(plotresultdir)) dir.create(plotresultdir)
    message(sprintf('Export plots of input imageObj for `%s` images.', 
                    paste0(export_plots, collapse = "`, `")))
  }
  
  if(length(morphToShow)!=3){
    stop("`morphToShow` must has 3 elements for channel index of R, B, G channels, respectively.")
  }
  
  if(!identical(dim(rgbMorphObj)[c(1,2,4)], dim(labelsObj))){
    stop("Mismatched dimenstion between `rgbMorphObj` and `labelsObj`.")
  }
  
  if(!all(morphToShow %in% seq_len(dim(rgbMorphObj)[3]))){
    stop("The input `morphToShow` is beyond the channel range of `rgbMorphObj`.")
  }
  
  
  if(!is.null(cell_ids)){
    if(length(cell_ids) != dim(labelsObj)[3]){
      stop("The provided `cell_ids` doesn't match with number of frame in input image obj.")
    }
    
    if(explotFlag){
      message("The provided `cell_ids` would be printed on top of exported plots.")
    }
  }
  
  # paint cell borders on RGB morphology 
  segImg <-  EBImage::Image(rgbMorphObj[, , morphToShow, ], colormode = EBImage::Color)
  segImg <-  EBImage::paintObjects(labelsObj, segImg, col= border_color)
  
  if(explotFlag){
    imagesPerRow <- round(sqrt(dim(labelsObj)[3]))
    
    # Get dimensions of combined image
    img.dim <- dim(labelsObj)
    
    if(!is.null(cell_ids)){
      xcenter <- img.dim[1]/2 + ((seq_len(img.dim[3])-1)%%imagesPerRow)*img.dim[1]
      yborder <- ceiling(seq_len(img.dim[3])/imagesPerRow-1)*img.dim[2]
      
    }
    
    
    if('labels' %in% export_plots){
      png(fs::path(plotresultdir, paste0(fileprefix, '__cellLabel.png')), 
          width = 3 * imagesPerRow, height = 3*ceiling(dim(labelsObj)[3]/imagesPerRow), 
          units = "in", res = 300)
      
      EBImage::display(EBImage::colorLabels(labelsObj), all = TRUE, 
                       method = 'raster', nx = imagesPerRow)
      if(!is.null(cell_ids)){
        rect(xcenter - img.dim[1]*0.4, yborder+img.dim[2]*0.2,
             xcenter + img.dim[1]*0.4, yborder,
             col = "#FFFFFF80", border = NA)
        text(x = xcenter, y = yborder, 
             labels = cell_ids, adj = c(1, 0.5), pos = 1, col = "black")
      }
      
      dev.off()
    }
    
    if('morph' %in% export_plots){
      png(fs::path(plotresultdir, paste0(fileprefix, '__morphRGB_', 
                                         paste0(morphToShow, collapse = "-"), 
                                         '.png')), 
          width = 3 * imagesPerRow, height = 3*ceiling(dim(labelsObj)[3]/imagesPerRow), 
          units = "in", res = 300)
      EBImage::colorMode(rgbMorphObj) <- EBImage::Color
      EBImage::display(rgbMorphObj[, , morphToShow, ], method = 'raster', 
                       all = TRUE, nx = imagesPerRow)
      if(!is.null(cell_ids)){
        rect(xcenter - img.dim[1]*0.4, yborder+img.dim[2]*0.2,
             xcenter + img.dim[1]*0.4, yborder,
             col = "#FFFFFF80", border = NA)
        text(x = xcenter, y = yborder, 
             labels = cell_ids, adj = c(1, 0.5), pos = 1, col = "black")
      }
      
      dev.off()
      
    }
    
    if('overlay' %in% export_plots){
      png(fs::path(plotresultdir, paste0(fileprefix, '__border_on_morph.png')), 
          width = 3 * imagesPerRow, height = 3*ceiling(dim(labelsObj)[3]/imagesPerRow), 
          units = "in", res = 300)
      
      EBImage::display(segImg, all = TRUE, 
                       method = 'raster', nx = imagesPerRow)
      if(!is.null(cell_ids)){
        rect(xcenter - img.dim[1]*0.4, yborder+img.dim[2]*0.2,
             xcenter + img.dim[1]*0.4, yborder,
             col = "#FFFFFF80", border = NA)
        text(x = xcenter, y = yborder, 
             labels = cell_ids, adj = c(1, 0.5), pos = 1, col = "black")
      }
      
      dev.off()
    }
  }
  
  return(segImg)
  
}


#' @title colorCellsByMeta
#' @description color cell labels with cell metadata 
#' @param labelsImgs named list of cell label image arrays with one query cell per element, each label image array is a 2D matrix with value in integer indicating cell id
#' @param fullCellMeta data.frame for cell metadata of all cells in study
#' @param cellID_coln column name of cell_ID in `fullCellMeta`
#' @param cellLabel_coln column name of cell label values in `fullCellMeta`
#' @param UID_coln column name of unique identifier for each fov in `fullCellMeta`
#' @param color_coln column name of metadata used to color cell labels in `fullCellMeta`
#' @param isDiscrete flag discrete metadata (default = TRUE)
#' @param colors_to_use a vector of colors to use; when isDiscrete = TRUE, must have length and names matched to unique values in color_coln 
#' @param value_limits a vector of two elements for the value range to display, used when isDiscrete = FALSE; default = NULL to use the 1% and 99% quantile of values in color_coln across entire study.
#' @param plotresultdir output folder for color legend used 
#' @return a list of single-frame rgb image array with one element per query cell
#' @export 
colorCellsByMeta <- function(labelsImgs, 
                             fullCellMeta, 
                             cellID_coln = 'cell_ID', 
                             cellLabel_coln = 'CellId',
                             UID_coln = 'UID',
                             color_coln = 'celltype',
                             isDiscrete = TRUE, 
                             colors_to_use = NULL, 
                             value_limits = NULL,
                             plotresultdir = getwd()){
  colns <- c(cellID_coln, cellLabel_coln, UID_coln, color_coln)
  if(!all(colns %in% colnames(fullCellMeta))){
    stop(sprintf("Missing essential columns in `fullCellMeta`: `%s`.", 
                 paste0(setdiff(colns, colnames(fullCellMeta)), collapse = "`, `")))
  }
  
  if(is.null(names(labelsImgs))){
    stop('The names of `labelsImgs` should be cell_ID in `fullCellMeta`.')
  }
  
  fullCellMeta[['queryCells']] <- (fullCellMeta[[cellID_coln]] %in% names(labelsImgs))
  
  fovs <- fullCellMeta[fullCellMeta[['queryCells']], UID_coln]
  if(length(fovs) != length(labelsImgs)){
    message(sprintf("`labelsImgs` has %d cells, but %d cells are not present in `fullCellMeta`, continue to process the remaining %d cells.", 
                    length(labelsImgs), length(labelsImgs) - length(fovs), length(fovs)))
    if(length(fovs)<1){
      stop("No cells shared bewteen `labelsImgs` and `fullCellMeta`.")
    }
  }
  
  
  plotDF <- fullCellMeta[fullCellMeta[[UID_coln]] %in% unique(fovs), c(colns, 'queryCells')]
  
  if(!dir.exists(plotresultdir)) dir.create(plotresultdir)
  
  if(isDiscrete){
    # paint cell label based on categorical value
    allCTs <- levels(as.factor(fullCellMeta[[color_coln]]))
    
    if(is.null(colors_to_use)){
      # first color is too dark in black background
      colors_to_use <- setNames(Seurat::DiscretePalette(length(allCTs)+1, 
                                                        palette = "polychrome")[-1], 
                                nm = allCTs)
    }else{
      if(is.null(names(colors_to_use))){
        if(length(colors_to_use) < length(allCTs)){
          stop(sprintf("%d colors provided but %d unique values in `color_coln` = %s.",
                       length(colors_to_use), length(allCTs), color_coln))
        }else{
          colors_to_use <- sample(colors_to_use, length(allCTs))
          colors_to_use <- setNames(colors_to_use, nm = allCTs)
        }
      }else{
        if(length(setdiff(allCTs, names(colors_to_use)))>1){
          stop(sprintf("Missing colors for unique values in `color_coln` = %s: `%s`.",
                       color_coln, paste0(setdiff(allCTs, names(colors_to_use)), collapse = "`, `")))
        }
        
        colors_to_use <- colors_to_use[allCTs]
      }
    }
    
    # convert the label to numeric given the order in color vector
    mypalette  <- setNames(seq_along(colors_to_use), 
                           nm = names(colors_to_use))
    valToPlot <- mypalette[plotDF[[color_coln]]]
    # assign any NA as 0, same as extracellular space 
    valToPlot[is.na(valToPlot)] <- 0
    
    # normalize to 0-1 for plotting 
    valToPlot <- valToPlot/length(colors_to_use)
    plotDF[['valToPlot']] <- valToPlot
    
    # set zero value to be black
    mypalette <- c("#000000", unname(colors_to_use))
    
    # plot color bar legend 
    png(fs::path(plotresultdir, paste0(color_coln, "_color_legend.png")), 
        res = 500, units = "in", height = 12, width = 6)
    par(mar = c(2, 2, 2, 5))
    plot(1:10, 1:10, type = "n", axes = FALSE, xlab = "", ylab = "", main = color_coln)
    legend("right", col = mypalette, 
           legend = c('extracellular', names(colors_to_use)), pch = 16, pt.cex = 2)
    dev.off()
    par(mar = c(5.1, 4.1, 4.1, 2.1))
    
  }else{
    # paint cell label based on numeric value
    
    # create new color palette, set zero value to be black
    if(is.null(colors_to_use)){
      colors_to_use <- c(rev(rainbow(10, start=rgb2hsv(col2rgb('cyan'))[1], 
                                     end=rgb2hsv(col2rgb('blue'))[1])), 
                         rev(rainbow(10, start=rgb2hsv(col2rgb('red'))[1], 
                                     end=rgb2hsv(col2rgb('yellow'))[1])))
    }
    mypalette <- colorRampPalette(colors_to_use)(2**8)
    mypalette <- c("#000000", mypalette)
    
    
    # normalize numeric value to 0~1 range, use 0 for no-cell range
    # display for 1~99% quantile of entire study
    if(is.null(value_limits)){
      value_limits <- quantile(fullCellMeta[[color_coln]], prob = c(0.01, 0.99))
    }else{
      value_limits <- sort(value_limits)[1:2]
    }
    
    valToPlot <- plotDF[[color_coln]]
    valToPlot <- (valToPlot - value_limits[1])/diff(value_limits)+0.01
    valToPlot <- pmax(0.01, pmin(valToPlot, 1))
    plotDF[['valToPlot']] <- valToPlot
    
    # plot color bar legend 
    png(fs::path(plotresultdir, paste0(color_coln, "_colorbar_legend.png")), 
        res = 500, units = "in", height = 12, width = 6)
    par(mar = c(0, 0, 3, 0))
    plot(0, xlim = c(0, 6), ylim = c(-0.5, 1.2), axes = FALSE, type = 'n',  
         main = color_coln)
    corrplot::colorlegend(mypalette[2:length(mypalette)],
                          labels = seq(value_limits[1], value_limits[2], 
                                       by = diff(value_limits)/9), 
                          vertical = TRUE,
                          xlim = c(2, 4), align = 'l', offset = 0)
    
    dev.off()
    par(mar = c(5.1, 4.1, 4.1, 2.1))
    
  }
  
  tmpDF <- split(plotDF, as.factor(plotDF[[UID_coln]]))
  imgData <- lapply(
    tmpDF, 
    function(df){
      cells_to_plot <- df[df[['queryCells']], cellID_coln]
      
      numConverter <- setNames(df[['valToPlot']], 
                               nm = as.character(df[[cellLabel_coln]]))
      # add 0 for extra cellular space 
      numConverter <- c("0" = 0, numConverter)
      
      # need to get value for whole affected fov 
      metaImgs <- lapply(
        cells_to_plot, 
        function(cellName){
          ll <- labelsImgs[[cellName]]
          # convert label to values
          ll <- array(numConverter[as.character(ll)], dim = dim(ll))
          
          # for any cells not present in metadata (likely due to not passing QC), plot same as extracellular space 
          ll[is.na(ll)] <- 0
          
          # covert values to rgb given colormap 
          ll <- EBImage::colormap(EBImage::Image(ll), palette = mypalette)
          return(ll)
        })
      names(metaImgs) <- cells_to_plot
      
      return(metaImgs)
    }
  )
  
  # get cell_id
  imgListCellIds <- unlist(lapply(imgData, function(x) names(x)))
  
  imgData <- do.call(c, imgData)
  names(imgData) <- imgListCellIds
  
  return(imgData)
  
}


#' @title combineMetaToMorph_imgList
#' @description combine 2 set of images into one set
#' @param rgbMorphImgs list of rgb image arrays with one query cell per element, each rgb image array has at least 3 channels with value in range of 0-1 for RGB color
#' @param rgbMetaImgs list of rgb image array with one query cell per element, each rgb image array has exactly 3 channels with zero value for extracellular space 
#' @param morphToShow a vector of 3 elements for the channel index of each RGB channel in `rgbMorphImgs`.
#' @return a multi-frame rgb image array with intracellular region colored by `rgbMetaImgs` but extracellular region colored by `rgbMorphImgs`
#' @export 
combineMetaToMorph_imgList <- function(rgbMorphImgs, 
                                       rgbMetaImgs, 
                                       morphToShow = c(1,2,3)){
  
  if(length(rgbMorphImgs) != length(rgbMetaImgs)){
    stop("Must have same length in `rgbMorphImgs`, `rgbMetaImgs`, with each element for one query cell.")
  }
  
  cells_to_use <- names(rgbMorphImgs)
  if(is.null(cells_to_use)){
    cells_to_use <- names(rgbMetaImgs)
    if(is.null(cells_to_use)){
      stop("No cell_ID provided in the names of either `rgbMorphImgs` or `rgbMetaImgs`.")
    }else{
      names(rgbMorphImgs) <- cells_to_use
    }
  }else{
    if(is.null(names(rgbMetaImgs))){
      names(rgbMetaImgs) <- cells_to_use
      message("No names in provided `rgbMetaImgs`, use the names of `rgbMorphImgs` as cell_ID.")
    }else{
      if(!all(names(rgbMetaImgs) %in% cells_to_use)){
        stop("The names of `rgbMorphImgs` and `rgbMetaImgs` are not matched.")
      }
    }
  }
  
  if(length(morphToShow)!=3){
    stop("`morphToShow` must has 3 elements for channel index of R, B, G channels, respectively.")
  }
  
  imgData <- lapply(
    cells_to_use, 
    function(cellName){
      # remove intracellular region in morph 
      combImgs <- rgbMorphImgs[[cellName]]
      
      if(!all(morphToShow %in% seq_len(dim(combImgs)[3]))){
        stop("The input `morphToShow` is beyond the channel range of `rgbMorphImgs`.")
      }
      combImgs <- combImgs[, , morphToShow]
      combImgs[rgbMetaImgs[[cellName]] !=0] <- 0
      
      # add meta info image to the intracellular region
      combImgs <- rgbMetaImgs[[cellName]] + combImgs
      combImgs <- EBImage::Image(combImgs)
      
      EBImage::colorMode(combImgs) <- EBImage::Color
      
      return(combImgs)
    })
  names(imgData) <- cells_to_use
  
  return(imgData)
  
}

#' @title displayGallery_imgObj
#' @description display multi-frame image object as cell gallery in plot window 
#' @param imgObj multi-frame image object for 2D cell label matrix of 3D rgb array 
#' @param cell_ids a vector of cell_id for each frame in input image obj (default = NULL if not to display cell_ids on plot) 
#' @param imgType type of input image in `imgObj`, labels for 2D matrix, morph for multi-channel 3D array
#' @param morphToShow a vector of 3 elements for the channel index of each RGB channel in `imgObj`.
#' @param method display method, 'raster' as gallery in plot, 'browser' as one cell per frame in viewer
#' @param allFrame flag to show all frames, refer to \code{EBImage::display}
#' @param colormode color mode to display, default = color
#' @param imagesPerRow number of images per row, default = NULL to pick automatically
#' @export 
displayGallery_imgObj <- function(imgObj, 
                                  cell_ids, 
                                  imgType = c('labels', 'morph'),
                                  morphToShow = c(1,2,3),
                                  method = c('raster', 'browser'), 
                                  allFrame = TRUE,
                                  colormode = c('color', 'grayscale'),
                                  imagesPerRow = NULL){
  imgType <- match.arg(imgType, c('labels', 'morph'))
  method <- match.arg(method, c('raster', 'browser'))
  colormode <- match.arg(colormode,  c('color', 'grayscale'))
  
  # in case of method = raster, allFrame = T, colormode = color, only 3 morph channels allowed 
  if(method == 'raster' && allFrame && colormode == 'color'){
    if(length(morphToShow)!=3){
      stop("`morphToShow` must has 3 elements for channel index of R, B, G channels, respectively.")
    }
  }
  
  
  
  # Get dimensions of combined image
  img.dim <- dim(imgObj)
  frameNum <- img.dim[length(img.dim)]
  
  if(is.null(imagesPerRow)){
    imagesPerRow <- round(sqrt(frameNum))
  }
  
  
  if(!is.null(cell_ids)){
    if(length(cell_ids) != frameNum){
      stop(sprintf("%d cells in `cell_ids`, notmatching %d frames in `imgObj`.", 
                   length(cell_ids), frameNum))
    }
  }
  
  imgObj <- EBImage::Image(imgObj, colormode = EBImage::Grayscale)
  if(colormode =='color'){
    EBImage::colorMode(imgObj) <- EBImage::Color
  }
  
  if(method == 'raster'){
    if(imgType == 'labels'){
      EBImage::display(EBImage::colorLabels(imgObj), all = allFrame, 
                       method = method, nx = imagesPerRow)
    }
    
    if(imgType == 'morph'){
      if(!all(morphToShow %in% seq_len(img.dim[3]))){
        stop("The input `morphToShow` is beyond the channel range of `imgObj`.")
      }
      EBImage::display(imgObj[, , morphToShow, ], method = method, 
                       all = allFrame, nx = imagesPerRow)
    }
    
    if(!is.null(cell_ids) && method == 'raster'){
      xcenter <- img.dim[1]/2 + ((seq_len(frameNum)-1)%%imagesPerRow)*img.dim[1]
      yborder <- ceiling(seq_len(frameNum)/imagesPerRow-1)*img.dim[2]
      
      rect(xcenter - img.dim[1]*0.4, yborder+img.dim[2]*0.2,
           xcenter + img.dim[1]*0.4, yborder,
           col = "#FFFFFF80", border = NA)
      text(x = xcenter, y = yborder, 
           labels = cell_ids, adj = c(1, 0.5), pos = 1, col = "black")
    }
    
    return(imgType)
    
  }
  
  # in case of viewer, must run the display line in the end and return nothing 
  if(method == 'browser'){
    if(imgType == 'labels'){
      EBImage::display(EBImage::colorLabels(imgObj), all = allFrame, 
                       method = method)
    }
    
    if(imgType == 'morph'){
      if(!all(morphToShow %in% seq_len(img.dim[3]))){
        stop("The input `morphToShow` is beyond the channel range of `imgObj`.")
      }
      
      EBImage::display(imgObj[, , morphToShow, ], method = method, 
                       all = allFrame)
      
    }
  }
}

#' @title displayCellId_singleImg
#' @description display cell_id of one single-frame label image, the identity of id for each label is shown as color in legend
#' @param imgLabel single-frame label image for one query cell, in format of a 2D matrix with value in integer indicating cell id
#' @param query_cellID cell_ID of input `labelsImg`
#' @param cellID_converter a named vector of how values in imgLabel is linked to cell_ID, useful when input is based on stitched image; default = NULL to show the label value as cell_ID
#' @param seed seed for randomized colors of cell_ID
#' @return a named vector of colors used for plotting cell labels with names for cell id in integer value 
#' @export 
displayCellId_singleImg <- function(imgLabel, 
                                    query_cellID = NULL,
                                    cellID_converter = NULL, 
                                    seed = 123){
  
  # get randomized colors 
  set.seed(seed)
  cells_to_plot <- sort(setdiff(unique(imgLabel), 0))
  
  colors_to_use <- replicate(3, 
                             sample(255, 
                                    length(cells_to_plot), 
                                    replace = ifelse(length(cells_to_plot)>255, TRUE, FALSE)))
  colors_to_use <- rgb(red = colors_to_use[, 1], 
                       green = colors_to_use[, 2], 
                       blue = colors_to_use[, 3], 
                       alpha = 255, maxColorValue = 255)
  
  # set zero value to be black
  mypalette <- c("#000000", colors_to_use)
  
  # convert the label to numeric given the order in color vector
  numConverter <- setNames(seq_along(colors_to_use)/length(colors_to_use), 
                           nm = as.character(cells_to_plot))
  
  # convert label to values
  cll <- array(numConverter[as.character(imgLabel)], dim = dim(imgLabel))
  
  # plot extracellular space as black 
  cll[is.na(cll)] <- 0
  
  # covert values to rgb given colormap 
  cll <- EBImage::colormap(EBImage::Image(cll), palette = mypalette)
  
  if(!is.null(cellID_converter)){
    # convert label value to name
    cells_to_plot2 <- cellID_converter[as.character(cells_to_plot)]
    if(any(is.na(cells_to_plot2))){
      message(sprintf("`imgLabel` has value not matching the names in the provided `cellID_converter`: `%s`. \nDisplay the naive label values and decoded ID in legend instead.",
                      paste0(cells_to_plot[is.na(cells_to_plot2)], collapse = "`, `")))
      tmpMeta <- decodeStitchedLabelVal(cells_to_plot[is.na(cells_to_plot2)])
      tmpMeta[['cell_ID']] <- paste0(tmpMeta[['fov']], '_', tmpMeta[['CellId']])
      
      
      cells_to_plot2[is.na(cells_to_plot2)] <- paste0(tmpMeta[['labVal']], " (", tmpMeta[['cell_ID']], ")")
    }
    
    cells_to_plot <- cells_to_plot2
  }
  
  
  # display cell label on left and then legend on right
  par(mfrow = c(1,2))
  par(mar = c(0, 0, 0, 0))
  EBImage::display(cll, method = 'raster')
  plot(0, 0, type = "n", axes = FALSE, xlab = "", ylab = "")
  legend("left", title = ifelse(is.null(query_cellID), "", query_cellID), 
         col = mypalette, legend = c('extracellular', cells_to_plot), pch = 16, pt.cex = 2)
  
  par(mfrow = c(1,1))
  par(mar = c(5.1, 4.1, 4.1, 2.1))
  
  return(setNames(unname(colors_to_use), 
                  nm = as.character(cells_to_plot)))
  
}
