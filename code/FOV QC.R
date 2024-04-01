## needed:
# - plots: functions called on "res"
# - flagging rule - own function, but called by runFOVQC



#' FOV QC 
#' 
#' Run the full FOV QC workflow:
#' 
#' @param counts Counts matrix
#' @param xy 2-column matrix of cells' xy positiong
#' @param fov Vector of FOV IDs
#' @param genes Vector of gene names (omit negprobes and falsecodes)
#' @param barcodes Vector of gene barcodes, aligned to genes
runFOVQC <- function(counts, xy, fov, genes, barcodes) {
  
  ## create a matrix of barcode bit expression over sub-FOV grids:
  # define grids, get per-square gene expression:
  gridinfo <- makeGrid(mat = counts, xy = xy, fov = fov, squares_per_fov = 49, min_cells_per_square = 25) 
  # convert to per-square barcode bit expression:
  bitcounts <- cellxgene2squarexbit(counts = counts, grid = gridinfo$gridid, genes = genes, barcodes = barcodes) 
  
  # normalize:
  bitcounts <- sweep(bitcounts, 1, rowSums(bitcounts), "/") * mean(rowSums(bitcounts))
  
  ## for every grid square, match it to "control" squares from other FOVs, and get its residuals from them:
  # get neighbors:
  comparators <- getNearestNeighborsByFOV(x = bitcounts, 
                                          fov = gridinfo$gridfov, 
                                          n_neighbors = 10)
  
  # get expected:
  yhat <- t(sapply(1:nrow(bitcounts), function(i) {
    colMeans(bitcounts[comparators[i, ], ])
  }))
  # get resids:
  resid = log2((bitcounts + 5) / (yhat + 5))
  rownames(resid) = rownames(bitcounts)
  
  ## summarize bias per FOV:
  fovstats <- summarizeFOVBias(resid, gridinfo$gridfov)
  
  return(list(fovstats = fovstats, resid = resid, gridinfo = gridinfo))
}


# plots
if (FALSE) {
  pheatmap(resid[order(gridinfo$gridfov[rownames(resid)]), ], cluster_rows = F,
           col = colorRampPalette(c("darkblue", "blue", "white", "red", "darkred"))(101),
           breaks = seq(-1,1,length.out = 100))
  
  bit = "c12B"
  plotinds <- !is.na(gridinfo$gridid)
  # plot expression:
  plot(xy, cex = 0.2, asp = 1,
       col = viridis_pal(option = "B")(121)[1 + pmin(bitcounts[match(gridinfo$gridid, rownames(resid)), bit], 120)], main = bit)
  # plot resids:
  plot(xy, cex = 0.2, asp = 1,
       col = colorRampPalette(c("darkblue", "blue", "grey80", "red", "darkred"))(101)[
         pmax(pmin(51 + resid[match(gridinfo$gridid, rownames(resid)), bit] * 50, 101), 1)], main = bit)
  
  pheatmap(fovstats$bias * (fovstats$p < 0.01),
           col = colorRampPalette(c("darkblue", "blue", "white", "red", "darkred"))(101),
           breaks = seq(-1,1,length.out = 100))
}


#' Get nearest neighbors, sampling diffusely across other FOVs. 
#' 
#' @param x Data matrix, obs in rows, variables in columns
#' @param fov vector of FOV IDs, aligned to rows of x
#' @param n_neighbors How many neighbors to record
#' @return A matrix of dimension nrow(obs) * n_neighbors, giving the row IDs for each obs's selected neighbors.
getNearestNeighborsByFOV <- function(x, fov, n_neighbors = 10) {
  # get nearest:
  topneighbors <- FNN::get.knnx(data = x, 
                                query = x, 
                                k = n_neighbors * 50)$nn.index
  # reduce neighbors: none from own FOV, max of max_per_fov from each other fov:
  reducedneighbors <- t(sapply(1:nrow(topneighbors), function(i) {
    tempneighbors <- topneighbors[i, ]
    # only accept 1 per FOV
    tempneighbors[duplicated(fov[topneighbors[i, ]])] <- NA
    # none from own FOV:
    tempneighbors[fov[topneighbors[i, ]] == fov[i]] <- NA
    return(tempneighbors[!is.na(tempneighbors)][1:n_neighbors])
  }))
  
  return(reducedneighbors)
}


#' Get gridded expression within FOVs
#' 
#' @param mat Expression matrix, linear-scale
#' @param xy 2-column matrix of xy locations
#' @param fov Vector of FOV IDs
#' @return A list: a matrix of expression in grid squares, a vector giving the grid square assignment of each cell.   
makeGrid <- function(mat, xy, fov, squares_per_fov = 49, min_cells_per_square = 25) {
  
  # make grids:
  ncuts <- floor(sqrt(squares_per_fov))
  grid <- rep(NA, nrow(mat))
  gridfov <- c()
  
  for (fovid in unique(fov)) {   
    inds <- fov == fovid
    grid[inds] <- paste0(fov[inds], "_", cut(xy[inds, 1], ncuts), "_", cut(xy[inds, 2], ncuts))
    # associate these grid IDs with this FOV:
    gridfov[unique(grid[inds])] <- rep(fovid, length(unique(grid[inds])))
  }
  
  # clean up low-data grids:
  toosmallgridids <- names(which(table(grid) < min_cells_per_square))
  grid[is.element(grid, toosmallgridids)] <- NA
  
  return(list(gridid = grid, gridfov = gridfov))
}


#' Convery cell x gene matrix to grid square * bit matrix
#' @param counts Counts matrix
#' @param grid Vector assigning cells to squares, aligned to rows of counts
#' @param genes Vector of gene names (omit negprobes and falsecodes)
#' @param barcodes Vector of gene barcodes, aligned to genes
cellxgene2squarexbit <- function(counts, grid, genes, barcodes) {
  
  # number of bits:
  nreportercycles <- nchar(barcodes[1]) / 2
  nbits <- nreportercycles * 4
  # parse barcodes:
  bitmat = matrix(0, length(setdiff(unique(grid), NA)), nbits)
  colnames(bitmat) <- paste0("c", rep(seq_len(nreportercycles), each = 4), c("B", "G", 'Y', "R"))
  rownames(bitmat) <- setdiff(unique(grid), NA)
  for (i in seq_len(nreportercycles)) {
    barcodeposition <- i*2
    barcodehere <- substr(barcodes, barcodeposition, barcodeposition)
    for (col in c("B", "Y", "G", "R")) {
      tempgenes <- setdiff(genes[barcodehere == col], NA)
      temptotal <- Matrix::rowSums(counts[, tempgenes])
      tempsquaretotal <- by(temptotal, grid, mean)
      bitmat[names(tempsquaretotal), paste0("c", i, col)] <- tempsquaretotal
    }
  }
  return(bitmat)
}

#' Convert gene counts to bit counts
#' @param mat Counts matrix
#' @param genes Vector of gene names (omit negprobes and falsecodes)
#' @param barcodes Vector of gene barcodes, aligned to genes
genes2bits <- function(mat, genes, barcodes) {
  
  # number of bits:
  nreportercycles <- nchar(barcodes[1]) / 2
  nbits <- nreportercycles * 4
  # parse barcodes:
  bitmat = matrix(0, nrow(mat), nbits)
  colnames(bitmat) = paste0("c", rep(seq_len(nreportercycles), each = 4), c("B", "G", 'Y', "R"))
  for (i in seq_len(nreportercycles)) {
    barcodeposition <- i*2
    barcodehere <- substr(barcodes, barcodeposition, barcodeposition)
    for (col in c("B", "Y", "G", "R")) {
      tempgenes <- genes[barcodehere == col]
      bitmat[, paste0("c", i, col)] = Matrix::rowSums(mat[, is.element(colnames(mat), tempgenes)])
    }
  }
  rownames(bitmat) <- rownames(mat)
  return(bitmat)
}

#' Summarize bias in FOVs
#' 
#' For each bit, and each FOV, get the mean bias and a p-value
summarizeFOVBias <- function(resid, gridfov) {
  gridfov = gridfov[rownames(resid)]
  fovs = unique(gridfov)
  bias <- p <- matrix(NA, length(fovs), ncol(resid),
                      dimnames = list(fovs, colnames(resid)))
  
  for (bit in colnames(resid)) {
    mod = summary(lm(resid[, bit] ~ as.factor(gridfov) - 1))$coef
    bias[gsub("as.factor\\(gridfov\\)", "", rownames(mod)), bit] = mod[, "Estimate"]
    p[gsub("as.factor\\(gridfov\\)", "", rownames(mod)), bit] = mod[, "Pr(>|t|)"]
  }
  
  return(list(bias = bias, p = p))
}

