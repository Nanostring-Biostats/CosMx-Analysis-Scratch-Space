#' Read flat files, return a sparse matrix of counts and a metadata data table
#' @param myflatfiledir Parent directory for flat files. Should hold one subdirectory of results for each slide;
#'   each of those subdirectories should have a exprMat_file.csv (or .csv.gz) and a metadata_file.csv (or .csv.gz).
#'   Of note:
#'   \itemize{
#'    \item{The possibly non-unique metadata column "fov" is replaced by the column "FOV", for which a given FOV ID will only appear in a single slide.}
#'    \item{The possibly non-unique metadata column "cell_ID" is deleted, leaving the column "cell_id", for which a given cell ID will only appear in a single slide & FOV.}
#'   }
#' @param slidenames Optional vector of which slide folders to read. If left NULL, all slide folders will be read.
#' @param return_negprobes Logical for whether to return the expression matrix of negprobes
#' @param return_negprobes Logical for whether to return the expression matrix of falsecodes
#' @param offset_slides Logical, for whether to shift slide positions so they're non-overlapping.
#' @param thisinstrument_nanometers_per_pixel Nanometers per pixel in your data. Default of 120.280945 applies to all commercial instruments. Your RunSummary file specifies the value for your instrument.
#' @return A list: "counts", "metadata", "xy", and optionally "negcounts" and "falsecounts" 
#' @importFrom data.table fread
#' @importFrom R.utils gunzip
#' @importFrom Matrix rowSums
readFlatFiles <- function(myflatfiledir, slidenames = NULL, return_negprobes = FALSE, return_falsecodes = FALSE, offset_slides = TRUE, thisinstrument_nanometers_per_pixel = 120.280945) {
  ### automatically get slide names:
  if (is.null(slidenames)) {
    slidenames <- dir(myflatfiledir)
  }
  
  ### lists to collect the counts matrices and metadata, one per slide
  countlist <- vector(mode='list', length=length(slidenames)) 
  metadatalist <- vector(mode='list', length=length(slidenames)) 
  
  for(i in 1:length(slidenames)){
    
    slidename <- slidenames[i] 
    
    msg <- paste0("Loading slide ", slidename, ", ", i, "/", length(slidenames), ".")
    message(msg)    
    # slide-specific files:
    thisslidesfiles <- dir(paste0(myflatfiledir, "/", slidename))
    
    # load in metadata:
    thisslidesmetadata <- thisslidesfiles[grepl("metadata\\_file.csv", thisslidesfiles)]
    if(length(thisslidesmetadata) != 1){stop(paste0("Check slide ", slidename, " for format issue. Not finding exactly one metadata file."))}
    tempdatatable <- data.table::fread(paste0(myflatfiledir, "/", slidename, "/", thisslidesmetadata))
    
    # numeric slide ID 
    slide_ID_numeric <- tempdatatable[1,]$slide_ID 
    
    # load in counts as a data table:
    thisslidescounts <- thisslidesfiles[grepl("exprMat\\_file", thisslidesfiles)]
    # if two files, take the one that's been unzipped:
    if (length(thisslidescounts) == 2) {
      thisslidescounts <- thisslidescounts[!grepl("\\.gz", thisslidescounts)]
    }
    if (length(thisslidescounts) == 0) {
      warning(paste0("No exprMat files found for ", slidename))
    }
    if (length(thisslidescounts) > 2) {
      warning(paste0("> 2 exprMat files found for ", slidename, '. Expecting 2 or less (zipped and unzipped)'))
    }
    # unzip counts file if needed:
    if (substr(thisslidescounts, nchar(thisslidescounts)-2, nchar(thisslidescounts)) == ".gz") {
      R.utils::gunzip(paste0(myflatfiledir, "/", slidename, "/", thisslidescounts), remove = FALSE, overwrite = TRUE)
      thisslidescounts <- gsub(".csv.gz", ".csv", thisslidescounts)
    }
    countsfile <- paste0(myflatfiledir, "/", slidename, "/", thisslidescounts)
    nonzero_elements_perchunk <- 5*10**7
    ### Safely read in the dense (0-filled ) counts matrices in chunks.
    lastchunk <- FALSE 
    skiprows <- 0
    chunkid <- 1
    
    required_cols <- fread(countsfile, select=c("fov", "cell_ID"))
    stopifnot("columns 'fov' and 'cell_ID' are required, but not found in the counts file" = 
                all(c("cell_ID", "fov") %in% colnames(required_cols)))
    number_of_cells <- nrow(required_cols)
    
    number_of_cols <-  ncol(fread(countsfile, nrows = 2))
    number_of_chunks <- ceiling(exp(log(number_of_cols) + log(number_of_cells) - log(nonzero_elements_perchunk)))
    chunk_size <- floor(number_of_cells / number_of_chunks)
    sub_counts_matrix <- vector(mode='list', length=number_of_chunks)
    
    pb <- txtProgressBar(min = 0, max = number_of_chunks, initial = 0, char = "=",
                         width = NA, title, label, style = 3, file = "")
    cellcount <- 0
    while(lastchunk==FALSE){
      read_header <- FALSE
      if(chunkid==1){
        read_header <- TRUE
      }
      
      countsdatatable <- data.table::fread(countsfile,
                                           nrows=chunk_size,
                                           skip=skiprows + (chunkid > 1),
                                           header=read_header)
      if(chunkid == 1){
        header <- colnames(countsdatatable)
      } else {
        colnames(countsdatatable) <- header
      }
      
      nr <- nrow(countsdatatable)
      cellcount <- cellcount + nr
      lastchunk <- (cellcount >= number_of_cells) || nr == 0L
      skiprows <- skiprows + nr
      
      # deleting next line: fov and id columns no longer in counts matrix:
      slide_fov_cell_counts <- paste0("c_",slide_ID_numeric, "_", countsdatatable$fov, "_", countsdatatable$cell_ID)
      sub_counts_matrix[[chunkid]] <- as(countsdatatable[,setdiff(colnames(countsdatatable), c("fov", "cell_ID")),with=FALSE], "sparseMatrix") 
      rownames(sub_counts_matrix[[chunkid]]) <- slide_fov_cell_counts 
      setTxtProgressBar(pb, chunkid)
      chunkid <- chunkid + 1
    }
    
    close(pb)   
    
    countlist[[i]] <- do.call(rbind, sub_counts_matrix) 
    # ensure that cell-order in counts matches cell-order in metadata   
    slide_fov_cell_metadata <- paste0("c_",slide_ID_numeric, "_", tempdatatable$fov, "_", tempdatatable$cell_ID)
    countlist[[i]] <- countlist[[i]][match(slide_fov_cell_metadata, rownames(countlist[[i]])),] 
    metadatalist[[i]] <- tempdatatable 
    
    # track common genes and common metadata columns across slides
    if(i==1){
      sharedgenes <- colnames(countlist[[i]]) 
      sharedcolumns <- colnames(tempdatatable)
    } else {
      lostgenes <- setdiff(sharedgenes, colnames(countlist[[i]]))
      if (length(lostgenes) > 0) {
        warnings(paste0("Dropping genes not present in other slides: ", paste0(lostgenes, collapse = ", ")))
      }
      lostgenes2 <- setdiff(colnames(countlist[[i]]), sharedgenes)
      if (length(lostgenes2) > 0) {
        warnings(paste0("Dropping genes not present in other slides: ", paste0(lostgenes2, collapse = ", ")))
      }
      sharedgenes <- intersect(sharedgenes, colnames(countlist[[i]]))
      sharedcolumns <- intersect(sharedcolumns, colnames(tempdatatable))
    }
    
  }
  
  # reduce to shared metadata columns and shared genes
  for(i in 1:length(slidenames)){
    metadatalist[[i]] <- metadatalist[[i]][, ..sharedcolumns]
    countlist[[i]] <- countlist[[i]][, sharedgenes]
  }
  
  counts <- do.call(rbind, countlist)
  metadata <- rbindlist(metadatalist)
  
  # isolate negative control matrices:
  negcounts <- counts[, grepl("Negative", colnames(counts))]
  falsecounts <- counts[, grepl("SystemControl", colnames(counts))]
  
  # reduce counts matrix to only genes:
  counts <- counts[, !grepl("Negative", colnames(counts)) & !grepl("SystemControl", colnames(counts))]
  
  # checks
  stopifnot(identical(rownames(counts), metadata$cell_id))
  stopifnot(max(abs(Matrix::rowSums(counts) - metadata$nCount_RNA)) == 0)
  
  ## process metadata ------------------------------------------------
  # add to metadata: add a global non-slide-specific FOV ID:
  metadata$FOV <- paste0("s", metadata$slide_ID, "f", metadata$fov)
  metadata$fov <- NULL
  
  # remove cell_ID metadata column, which only identifies cell within slides, not across slides:
  metadata$cell_ID <- NULL
  
  # add mean negprobe to metadata:
  metadata$negmean <- Matrix::rowMeans(negcounts)
  
  # extract xy positions as a separate object:
  xy <- as.matrix(metadata[, c("CenterX_global_px", "CenterY_global_px")])
  rownames(xy) <- metadata$cell_id
  # rescale to mm:
  xy <- xy * thisinstrument_nanometers_per_pixel / 1000000
  colnames(xy) <- paste0(c("x", "y"), "_mm")
  
  # make slides' xy non-overlapping:
  if (offset_slides) {
    xy <- condenseTissues(xy = xy, 
                          tissue = metadata$Run_Tissue_name, 
                          tissueorder = NULL,  # optional, specify the order tissues are tiled in
                          buffer = 1, # space between tissues
                          widthheightratio = 1) # desired shape of final xy locations 
  }
  
  ## returns:
  out <- list(counts = counts, metadata = metadata, xy = xy)
  if (return_negprobes) {
    out$negcounts <- negcounts
  }
  if (return_falsecodes) {
    out$falsecounts <- falsecounts
  }
  return(out)
}




#' Condense tissues' xy positions
#' (has a number of features, but only used in this pipeline to shift whole slides' xy coords to avoid overlap.)
condenseTissues <- function(xy, tissue, tissueorder = NULL, buffer = 0.2, widthheightratio = 4/3) {
  
  # get each tissue's dimensions:
  tissdf <- data.frame(tissue = unique(tissue))
  tissdf$width <- sapply(unique(tissue), function(tiss) {
    diff(range(xy[tissue == tiss, 1]))
  })
  tissdf$height <- sapply(unique(tissue), function(tiss) {
    diff(range(xy[tissue == tiss, 2]))
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
  tissuesperrow <- round(sqrt(nrow(tissdf)) * widthheightratio * mean(tissdf$height) / mean(tissdf$width))
  targetwidth <- sum(tissdf$width[1:tissuesperrow], na.rm = TRUE) + buffer * (tissuesperrow - 1)
  
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
