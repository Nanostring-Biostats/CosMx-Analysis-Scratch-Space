
#' Derive a polygon from transcript locations
#' 
#' Input x and y positions; derive an encircling polygon
#' @param x Vector of x positions or 2-column matrix
#' @param y Vector of y positions
#' @param type Either "ahull" or "chull". The latter is faster but limited to convex shapes.
#' @param alpha Passed to alphahull::ahull(). Smaller makes for more nuanced polygons.
#' @return A polygon's coordinated (a 2-column matrix)
#' @importFrom alphahull ahull
#' @importFrom grDevices chull
getPolyFromTranscripts <- function(x, y = NULL, type = "ahull", alpha = 0.1) {
  
  if (type == "ahull") {
    # get the points of the alpha hull:
    temp <- alphahull::ahull(x = x, y = y, alpha = alpha)
    out <- temp$ashape.obj$edges
    
    # place them in order by their polar coords:
    center = c(median(range(out[, "x1"])), median(range(out[, "y1"])))
    newx <- out[, "x1"] - median(range(out[, "x1"]))
    newy <- out[, "y1"] - median(range(out[, "y1"]))
    theta <- atan2(newy, newx)
    
    # return a polygon:
    return(out[c(order(theta), order(theta)[1]), c("x1", "y1")])
  } else if (type == "chull") {
    out <- chull(x, y)
    out <- c(out, out[1])
    if (!is.null(dim(x))) {
      print("2")
      return(as.matrix(x[out, ]))
    } else {
      return(cbind(x, y)[out, ])
    }
  }
  
}



