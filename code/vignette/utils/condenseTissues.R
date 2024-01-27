#' Arrange tissues in xy space to reduce whitespace
#' 
#' Uses a shelf algorithm: places tallest tissues on the bottom shelf, and so on. 
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
  targetwidth <- sum(tissdf$width[1:tissuesperrow]) + buffer * (tissuesperrow - 1)
  
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
    print(tissdf)
  }
  # now update xy:
  for (tiss in unique(tissue)) {
    inds <- tissue == tiss
    xy[inds, 1] <- xy[inds, 1] - min(xy[inds, 1]) + tissdf$x[tissdf$tissue == tiss]
    xy[inds, 2] <- xy[inds, 2] - min(xy[inds, 2]) + tissdf$y[tissdf$tissue == tiss]
  }
  return(xy)  
}