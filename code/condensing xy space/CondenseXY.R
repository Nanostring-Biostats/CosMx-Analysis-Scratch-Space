
# demo and test:
if (FALSE) {
  library(RColorBrewer)
  tissuexy = cbind(c(1,1,2,
                     5,6,
                     8,8,8,
                     1,2,3,4,4,
                     6),
                   c(1,2,2,
                     3.5,3.5,
                     0,1,2,
                     -2,-2,-2,-2,-1,
                     -1))
  tissue = fov = xy = c()
  for (i in 1:nrow(tissuexy)) {
    xy = rbind(xy, 
               sweep(cbind(runif(1000), runif(1000)), 2, tissuexy[i, ], "+"))
    fov = c(fov, rep(i, 1000))
  }
  tissue[1:3000] = "a"
  tissue[3001:5000] = "a"
  tissue[5001:8000] = "a"
  tissue[8001:13000] = "b"
  tissue[13001:14000] = "b"
  
  
  cols = c(brewer.pal(8, "Set3"), brewer.pal(8, "Set2"))
  par(mfrow = c(3,2))
  par(mar = c(0,0,0,0))
  plot(xy, col = cols[as.numeric(as.factor(tissue))], asp = 1)
  plot(xy, col = cols[as.numeric(as.factor(fov))], asp = 1)
  
  
  newxy = condenseXY(xy = xy, fov = fov, tissue = tissue, 
                     condenseFOVs = TRUE, condensetissues = TRUE, # which operations to perform
                     eps = 1.5, mindist = 0.1, fovbuffer = 0.2, # FOV condensing args
                     tissueorder = NULL, tissuebuffer = 0.2, widthheightratio = 4/3) # Tissue condensing args
  par(mfrow = c(1,2))
  plot(xy, col = cols[as.numeric(as.factor(tissue))], asp = 1)
  plot(newxy, col = cols[as.numeric(as.factor(tissue))], asp = 1)
}


#### main xy space condensing function ---------------


#' Wrapper for within- and between-tissue space condensing
condenseXY <- function(xy, fov, tissue, 
                       condenseFOVs = TRUE, condensetissues = TRUE, # which operations to perform
                       eps = 1.5, mindist = 0.1, fovbuffer = 0.2, # FOV condensing args
                       tissueorder = NULL, tissuebuffer = 0.2, widthheightratio = 4/3) {  # Tissue condensing args
  if (condenseFOVs) {
    xy <- condenseFOVs(xy = xy, fov = fov, tissue = tissue, eps = eps, mindist = mindist, buffer = fovbuffer) 
  }
  if (condensetissues) {
    xy <- condenseTissues(xy = xy, tissue = tissue, tissueorder = tissueorder, buffer = tissuebuffer, widthheightratio = widthheightratio) 
  }
  return(xy)  
}


#### Functions for condensing FOVs --------------------------------

#' Condense FOVs from each of many tissues
condenseFOVs <- function(xy, fov, tissue, eps = 1.5, mindist = 0.1, buffer = 0.2) {
  
  newxylist <- sapply(unique(tissue), function(tiss) {
    inds <- tissue == tiss
    condenseFOVs_onetissue(xy = xy[inds, ], 
                           fov = fov[inds], 
                           eps = eps, 
                           mindist = mindist, 
                           buffer = buffer)
  })
  newxy <- xy * NA
  for (tiss in unique(tissue)) {
    inds <- tissue == tiss
    newxy[inds, ] <- newxylist[[tiss]]
  }
  return(newxy)
}



condenseFOVs_onetissue <- function(xy, fov, eps = 1.5, mindist = 0.1, buffer = 0.2, ntries = 3) {
  
  # define FOV groups:
  fovcenters <- t(sapply(unique(fov), function(fid) {
    inds <- fov == fid
    return(midpoint(xy[inds, ]))
  }))
  # give each fov a group:
  fovgroup <- dbscan::dbscan(fovcenters, eps = eps, minPts = 1)$cluster
  fovgroup[fovgroup == 0] <- max(fovgroup) + 1:sum(fovgroup == 0)
  # give each cell a group:
  group <- fovgroup[match(fov, unique(fov))]
  
  # get bounding points of each FOV group:
  isboundary <- TRUE
  
  # get order in which to move groups:
  center <- midpoint(xy)
  groupdists <- sapply(unique(group), function(groupid) {
    inds <- group == groupid
    tempcenter <- midpoint(xy[inds, ])
    return((sum(tempcenter - center)^2))
  })
  grouporder <- unique(group)[order(groupdists)]
  propmove <- 0.5
  
  # initialize results for random paths:
  xy0 <- xy
  solutions <- list()
  # try a few solution paths:
  for (try in 1:ntries) {
    xy <- xy0
    center <- midpoint(xy)
    
    # pull blocks towards the center:
    for (groupid in rep(grouporder, 3)) {
      
      thisgroup <- group == groupid
      direction <- sample(1:3, 1)
      if (direction == 1) {
        # move directly towards the center:
        vec <- center - midpoint(xy[thisgroup, ])
        theta <- atan2(y = vec[2], x = vec[1]) 
        movedist <- maxMove(xy = xy[thisgroup & isboundary, ], 
                            xy2 = rbind(center, xy[!thisgroup & isboundary, ]),
                            vec = vec, 
                            mindist = mindist) * propmove
        xy[thisgroup, 1] <- xy[thisgroup, 1] + cos(theta) * pmax(movedist - buffer, 0)
        xy[thisgroup, 2] <- xy[thisgroup, 2] + sin(theta) * pmax(movedist - buffer, 0)
      }
      if (direction == 2) {
        # move up towards the center:
        vec <- center - midpoint(xy[thisgroup, ])
        vec <- c(0, vec[2])
        theta <- atan2(y = vec[2], x = vec[1]) 
        movedist <- maxMove(xy = xy[thisgroup & isboundary, ], 
                            xy2 = rbind(center, xy[!thisgroup & isboundary, ]),
                            vec = vec, 
                            mindist = mindist) * propmove
        xy[thisgroup, 1] <- xy[thisgroup, 1] + cos(theta) * pmax(movedist - buffer, 0)
        xy[thisgroup, 2] <- xy[thisgroup, 2] + sin(theta) * pmax(movedist - buffer, 0)
      }
      if (direction == 3) {
        # move sideways towards the center:
        vec <- center - midpoint(xy[thisgroup, ])
        vec <- c(vec[1], 0)
        theta <- atan2(y = vec[2], x = vec[1]) 
        movedist <- maxMove(xy = xy[thisgroup & isboundary, ], 
                            xy2 = rbind(center, xy[!thisgroup & isboundary, ]),
                            vec = vec, 
                            mindist = mindist) * propmove
        xy[thisgroup, 1] <- xy[thisgroup, 1] + cos(theta) * pmax(movedist - buffer, 0)
        xy[thisgroup, 2] <- xy[thisgroup, 2] + sin(theta) * pmax(movedist - buffer, 0)
      }
      plot(xy, col = cols[fov], asp = 1)
      
      center <- midpoint(xy)
      
      if (groupid == grouporder[length(grouporder)]) {
        propmove = mean(c(propmove, 1))
      }
      
    }
    solutions[[try]] <- xy  
  }
  # choose best version:
  boundingAreas <- sapply(solutions, boundingArea)
  
  return(solutions[[which.min(boundingAreas)]])
}

midpoint <- function(xy) {
  apply(xy, 2, function(x){median(range(x))})
}


#' dist and direction from a set of points (presumably a bounding polygon) to a point xy0:
minDist <- function(xy, xy0) {
  dists <- sweep(xy, 2, xy0, "-")^2
  mindist <- min(dists)
  closestpoint <- which.min(dists)
  diffvec <- xy0 - xy[closestpoint, ]
  return(list(mindist = mindist, diffvec = diffvec))
}


rotateCoords <- function(xy, vec) {
  theta <- atan2(y = vec[2], x = vec[1]) 
  newxy <- xy * NA
  newxy[, 1] <- xy[, 1] * cos(theta) + xy[, 2] * sin(theta)
  newxy[, 2] <- - xy[, 1] * sin(theta) + xy[, 2] * cos(theta)
  return(newxy)
}

#' max dist before a point collides:
getMaxMove <- function(xy0, xy2, mindist) {
  # subset of points that are close to the vec projected from xy[i, ], and along it instead of behind it:
  toconsider <- (abs(xy2[, 2] - xy0[2]) < mindist) & (xy2[,1] > xy0[1])
  if (sum(toconsider) > 0) {
    return(min(xy2[toconsider, 1] - xy0[1]))
  } else {
    return(NA)
  }
} 


#' how far can a set of points move in a direction before hitting others?
maxMove <- function(xy, xy2, vec, mindist) {
  
  # transform coordinate space so vec is horizontal:
  xy <- rotateCoords(xy, vec)
  xy2 <- rotateCoords(xy2, vec)
  
  # for every point in xy, find how far it can go in vec until it hits something:
  dist2collison <- apply(xy, 1, getMaxMove, xy2, mindist)
  maxdist <- suppressWarnings(min(dist2collison, na.rm = T))
  if (maxdist == Inf) {
    maxdist <- 0
  }
  return(maxdist)
  
}




#' measure area of bounding box (used to choose the most condensed solution):
boundingArea <- function(xy) {
  diff(range(xy[, 1])) * diff(range(xy[, 2]))
}




#' Choose plot dimensions given a range in xy space:
#' Assumes points drawn with cex = 0.1
xy2inches <- function(span) {
  return(span * 0.7)  # well-considered scaling factor
}


#### Tissue condensing functions ---------------


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



getdimensions <- function(xy) {
  width = diff(range(xy[,1]))
  height = diff(range(xy[,2]))
}

