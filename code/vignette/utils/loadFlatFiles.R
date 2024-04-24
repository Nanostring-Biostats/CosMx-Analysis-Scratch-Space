
#' Load counts matrix from gzipped flat-file
#'
#' Safely read in the dense (0-filled ) counts matrices in chunks.
#' Expected format: countsfile should point to a cells x genes countsfile matrix including cell-identifier columns named "fov" and "cell_ID".  
#' The remaining columns should be one per gene, indicating the counts of that gene for each cell.
#' @param countsfile possibly gzipped countsfile csv.
#' @param slide_ID_numeric numeric slide ID used to create cell_ID.
#' @param nonzero_elements_perchunk  
#' @return A sparse cells x genes counts matrix.
read_gzipped_countsfile <- function(countsfile, slide_ID_numeric=1, nonzero_elements_perchunk=5*10**7){
    
  ### Safely read in the dense (0-filled ) counts matrices in chunks.
  ### Note: the file is gzip compressed, so we don't know a priori the number of chunks needed.
  lastchunk <- FALSE 
  skiprows <- 0
  chunkid <- 1
  
  required_cols <- fread(countsfile, select=c("fov", "cell_ID"))
  stopifnot("columns 'fov' and 'cell_ID' are required, but not found in the counts file" = 
              all(c("cell_ID", "fov") %in% colnames(required_cols)))
  number_of_cells <- nrow(required_cols)
  
  number_of_cols <-  ncol(fread(countsfile, nrows = 2))
  number_of_chunks <- ceiling(number_of_cols * number_of_cells / (nonzero_elements_perchunk))
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
    
    countsdatatable <- data.table::fread(countsfile
                                         ,nrows=chunk_size
                                         ,skip=skiprows + (chunkid > 1)
                                         ,header=read_header
                                         )
    if(chunkid == 1){
      header <- colnames(countsdatatable)
    } else {
      colnames(countsdatatable) <- header
    }
     
    cellcount <- nrow(countsdatatable) + cellcount     
    if(cellcount == number_of_cells) lastchunk <- TRUE
    skiprows <- skiprows + chunk_size
    slide_fov_cell_counts <- paste0("c_",slide_ID_numeric, "_", countsdatatable$fov, "_", countsdatatable$cell_ID)
    sub_counts_matrix[[chunkid]] <- as(countsdatatable[,-c("fov", "cell_ID"),with=FALSE], "sparseMatrix") 
    rownames(sub_counts_matrix[[chunkid]]) <- slide_fov_cell_counts 
    setTxtProgressBar(pb, chunkid)
    chunkid <- chunkid + 1
  }
  
  close(pb)   
  
  return(do.call(rbind, sub_counts_matrix))
}




#' Load flat files exported from AtoMx
#' 
#' Given a directory holding flat files exported by AtoMx, loads them in. 
#' If flat files from multiple slides are present, merges them, preserving all shared genes and metadata variables.
#' 
#' Expected structure for flat files: exactly as exported by AtoMx: an overarching folder 
#' holds a folder for each slide. Each slide's folder contains at a minimum: 
#' an "exprMat" .csv file holding raw counts and a "metadata" .csv file holding cell metadata. 
#' 
#' @param flatfiledir Directory where flat files are located. 
#' @param countsflatfile_nonzero_elements_perchunk maximum number of nonzero elements to read in at a time from large counts matrix flat file. 
#' @return A list with 4 elements: \code{counts}, a sparse matrix of cell x gene count values, 
#'  and \code{metadata}, a data table of cell metadata, \code{negcounts}, a sparse matrix of cell x negprobe count values, 
#'  \code{falsecounts} a sparse matrix of cell x falsecode count values.
loadAtoMxFlatFiles <- function(flatfiledir, countsflatfile_nonzero_elements_perchunk=5*10**7) {
  ### automatically get slide names:
  slidenames <- dir(flatfiledir)


  ### lists to collect the counts matrices and metadata, one per slide
  countlist <- vector(mode='list', length=length(slidenames)) 
  metadatalist <- vector(mode='list', length=length(slidenames)) 
  
  for(i in 1:length(slidenames)){
     
    slidename <- slidenames[i] 

    msg <- paste0("Loading slide ", slidename, ", ", i, "/", length(slidenames), ".")
    message(msg)    
    # slide-specific files:
    thisslidesfiles <- dir(paste0(flatfiledir, "/", slidename))
    
    # load in metadata:
    thisslidesmetadata <- thisslidesfiles[grepl("metadata\\_file", thisslidesfiles)]
    tempdatatable <- data.table::fread(paste0(flatfiledir, "/", slidename, "/", thisslidesmetadata))
    
    # numeric slide ID 
    slide_ID_numeric <- tempdatatable[1,]$slide_ID 
      
    # load in counts as a data table:
    thisslidescounts <- thisslidesfiles[grepl("exprMat\\_file", thisslidesfiles)]
    countlist[[i]] <- read_gzipped_countsfile(paste0(flatfiledir, "/", slidename, "/", thisslidescounts),
                                              slide_ID_numeric,
                                              countsflatfile_nonzero_elements_perchunk
                                              ) 
    
    # ensure that cell-order in counts matches cell-order in metadata   
    slide_fov_cell_metadata <- paste0("c_",slide_ID_numeric, "_", tempdatatable$fov, "_", tempdatatable$cell_ID)
    countlist[[i]] <- countlist[[i]][match(slide_fov_cell_metadata, rownames(countlist[[i]])),] 
    metadatalist[[i]] <- tempdatatable 
    
    ##  **sanity check , loading all cells at once version should match**  ## # numeric slide ID 
    ##  **sanity check , loading all cells at once version should match**  ## slide_ID_numeric <- tempdatatable[1,]$slide_ID 
    ##  **sanity check , loading all cells at once version should match**  ## 
    ##  **sanity check , loading all cells at once version should match**  ## # create cell ID 
    ##  **sanity check , loading all cells at once version should match**  ## slide_fov_cell_counts <- paste0("c_",slide_ID_numeric, "_", countsdatatable$fov, "_", countsdatatable$cell_ID)
    ##  **sanity check , loading all cells at once version should match**  ## counts_matrix <- as(countsdatatable[,-c("fov", "cell_ID"),with=FALSE], "sparseMatrix") 
    ##  **sanity check , loading all cells at once version should match**  ## rownames(counts_matrix) <- slide_fov_cell_counts
    
    # track common genes and common metadata columns across slides
    if(i==1){
      sharedgenes <- colnames(countlist[[i]]) 
      sharedcolumns <- colnames(tempdatatable)
    }  else {
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
  
  # add to metadata: replace slide-specific FOV ID with a unique FOV ID:
  metadata$globalFOV <- paste0("s", metadata$slide_ID, "f", metadata$fov)
  
  # remove cell_ID metadata column, which only identifies cell within slides, not across slides:
  metadata$cell_ID <- NULL
  
  # isolate negative control matrices:
  negcounts <- counts[, grepl("Negative", colnames(counts))]
  falsecounts <- counts[, grepl("SystemControl", colnames(counts))]
  
  # reduce counts matrix to only genes:
  counts <- counts[, !grepl("Negative", colnames(counts)) & !grepl("SystemControl", colnames(counts))]
  
  return(list(counts = counts, metadata = metadata, negcounts = negcounts, falsecounts = falsecounts))
}