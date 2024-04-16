#' Draw polygons
#' 
#' Draws polygons for given cell IDs, calculating them from the transcript data frame if needed
#' @param polys List of polygons
#' @param cell_ids Vector of cell IDs
#' @param col Vector of colors, aligned to the cell_IDs
#' @param border Vector of border colors, aligned to the cell_IDs
#' @param transcript_df Data frame of transcript info. Needs these columns: "cell_id", "x", "y".
#' @param type Either "ahull" or "chull". The latter is faster but limited to convex shapes.
#' @param alpha Passed to alphahull::ahull(). Smaller makes for more nuanced polygons.
#' @param add Logical, for whether to open a new plot. If FALSE then polys are drawn atop existing plot.
#' @param outputpolys Whether to return the polygons list (useful if you've calculated more polygons within the function call). 
#' @param ... Arguments passed to plot(), in the event that add == FALSE
#' @return A list of polygons with the entries for the cell IDs filled in
#' @export
#' @importFrom graphics polygon
drawPolys <- function(polys, cell_ids = NULL, col = "#00008B80", border = NA, 
                      transcript_df = NULL, type = "chull", alpha = 0.1,
                      add = FALSE, outputpolys = FALSE, ...) {
  # draw all cells if no subset specified:
  if (is.null(cell_ids)) {
    cell_ids <- names(polys)
  }
  # derive any missing polys:
  nullpolys <- which(sapply(polys, is.null))
  needpolys <- intersect(names(nullpolys), cell_ids)
  if (length(needpolys) > 0) {
    message(paste0("deriving ", length(needpolys), " cells\' polgons"))
  }
  polys <- cellPolys(polys = polys, 
                     cell_ids = needpolys, 
                     transcript_df = transcript_df, 
                     overwrite = FALSE, 
                     type = type, alpha = alpha)
  # start a new plot if needed:
  if (!add) {
    xlim <- range(unlist(sapply(polys, function(xy){xy[,1]})))
    ylim <- range(unlist(sapply(polys, function(xy){xy[,2]})))
    plot(c(0,0), col = 0, xlim = xlim, ylim = ylim, ...)
  }
  # draw each cell's poly:
  if (length(col) == 1) {
    col <- rep(col, length(cell_ids))
  }
  if (length(border) == 1) {
    border <- rep(border, length(cell_ids))
  }
  for (i in 1:length(cell_ids)) {
    graphics::polygon(polys[[cell_ids[i]]], col = col[i], border = border[i])
  }
  if (outputpolys) {
    return(polys)
  }
}


