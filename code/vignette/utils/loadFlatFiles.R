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
loadAtoMxFlatFiles <- function(dir) {
  ### automatically get slide names:
  slidenames <- dir(flatfiledir)

  #### load in metadata from each slide:
  metadatalist <- sapply(slidenames, function(slidename) {
    thisslidesfiles <- dir(paste0(flatfiledir, "/", slidename))
    thisslidesmetadata <- thisslidesfiles[grepl("metadata\\_file", thisslidesfiles)]
    tempdatatable <- data.table::fread(paste0(flatfiledir, "/", slidename, "/", thisslidesmetadata))
    tempdatatable$slidename <- slidename
    return(list(tempdatatable))
  })
  
  ### harmonize metadata column names:
  # get column names shared by all metadata files:
  sharedcolumns <- colnames(metadatalist[[1]])
  if (length(metadatalist) > 1) {
    for (i in 2:length(metadatalist)) {
      sharedcolumns <- intersect(sharedcolumns, colnames(metadatalist[[i]]))
    }
  }
  # reduce to only shared columns:
  for (i in 1:length(metadatalist)) {
    metadatalist[[i]] <- metadatalist[[i]][, ..sharedcolumns]
  }
  
  ### load in counts matrix from each slide:
  countlist <- sapply(slidenames, function(slidename) {
    thisslidesfiles <- dir(paste0(flatfiledir, "/", slidename))
    thisslidescounts <- thisslidesfiles[grepl("exprMat\\_file", thisslidesfiles)]
    # load in counts as a data table:
    countsdatatable <- data.table::fread(paste0(flatfiledir, "/", slidename, "/", thisslidescounts))
    return(list(countsdatatable))
  })
  
  ### get shared genes:
  sharedgenes <- colnames(countlist[[1]])
  if (length(countlist) > 1) {
    for (i in 2:length(countlist)) {
      sharedgenes <- intersect(sharedgenes, colnames(countlist[[i]]))
    }
  }
  sharedgenes <- setdiff(sharedgenes, c("fov", "cell_ID"))
  
  #### align counts matrices to match metadata cell IDs, and convert to sparse matrices:
  for (i in 1:length(countlist)) {
    # align:
    metadata_cell_fov <- paste0(metadatalist[[i]]$fov, "_", metadatalist[[i]]$cell_ID)
    counts_cell_fov <- paste0(countlist[[i]]$fov, "_", countlist[[i]]$cell_ID)
    countlist[[i]] <- countlist[[i]][match(metadata_cell_fov, counts_cell_fov), ]
    # trim redundant columns:
    countlist[[i]] <- countlist[[i]][, ..sharedgenes]
    # convert to sparse matrix:
    countlist[[i]] <- as(countlist[[i]], "sparseMatrix")
  }
  
  ### condense metadata to a single data table:
  metadata <- c()
  counts <- c()
  
  for (i in 1:length(countlist)) {
    counts <- rbind(counts, countlist[[i]])
    metadata <- rbind(metadata, metadatalist[[i]])
  }
  # add cell IDs to counts:
  if(any(duplicated(metadata$cell_id))) {
    stop("Found duplicated cell IDs, probably from different slides. Make sure they're all unique.")
  }
  rownames(counts) <- metadata$cell_id
  colnames(counts) <- sharedgenes
  
  # add to metadata: replace slide-specific FOV ID with a unique FOV ID:
  metadata$FOV <- paste0("s", as.numeric(as.factor(metadata$slidename)), "f", metadata$fov)
  metadata$fov <- NULL
  
  # remove cell_ID metadata column, which only identifies cell within slides, not across slides:
  metadata$cell_ID <- NULL
  
  # isolate negative control matrices:
  negcounts <- counts[, grepl("Negative", colnames(counts))]
  falsecounts <- counts[, grepl("SystemControl", colnames(counts))]
  # reduce counts matrix to only genes:
  counts <- counts[, !grepl("Negative", colnames(counts)) & !grepl("SystemControl", colnames(counts))]
  
  return(list(counts = counts, metadata = metadata, negcounts = negcounts, falsecounts = falsecounts))
}