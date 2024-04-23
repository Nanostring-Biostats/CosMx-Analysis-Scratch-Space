#' Load flat files exported from AtoMx
#' 
#' Given a directory holding flat files exported by AtoMx, loads them in. 
#' If flat files from multiple slides are present, merges them, preserving all shared genes and metadata variables.
#' 
#' Expected structure for flat files: exactly as exported by AtoMx: an overarching folder 
#' holds a folder for each slide. Each slide's folder contains at a minimum: 
#' an "exprMat" .csv file holding raw counts and a "metadata" .csv file holding cell metadata. 
#' 
#' @param dir Directory where flat files are located. 
#' @return A list with 2 elements: \code{counts}, a sparse matrix of cell x gene count values, 
#'  and \code{metadata}, a data table of cell metadata.
loadAtoMxFlatFiles <- function(flatfiledir) {
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
      
    # load in counts as a data table:
    thisslidescounts <- thisslidesfiles[grepl("exprMat\\_file", thisslidesfiles)]
    countsdatatable <- data.table::fread(paste0(flatfiledir, "/", slidename, "/", thisslidescounts))
   
    # numeric slide ID 
    slide_ID_numeric <- tempdatatable[1,]$slide_ID 
    
    # create cell ID 
    slide_fov_cell_counts <- paste0("c_",slide_ID_numeric, "_", countsdatatable$fov, "_", countsdatatable$cell_ID)
    counts_matrix <- as(countsdatatable[,-c("fov", "cell_ID"),with=FALSE], "sparseMatrix") 
    rownames(counts_matrix) <- slide_fov_cell_counts
 
    # ensure that cell-order in counts matches cell-order in metadata   
    slide_fov_cell_metadata <- paste0("c_",slide_ID_numeric, "_", tempdatatable$fov, "_", tempdatatable$cell_ID)
    counts_matrix <- counts_matrix[match(slide_fov_cell_metadata, rownames(counts_matrix)),] 
   
    metadatalist[[i]] <- tempdatatable 
    countlist[[i]] <- counts_matrix 
    
    # track common genes and common metadata columns across slides
    if(i==1){
      sharedgenes <- colnames(counts_matrix) 
      sharedcolumns <- colnames(tempdatatable)
    }  else {
      sharedgenes <- intersect(sharedgenes, colnames(counts_matrix))
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