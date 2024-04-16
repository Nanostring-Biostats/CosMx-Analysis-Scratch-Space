
#' Initialize an empty polygons list
#' 
#' @param cell_ids Vector of cell IDs
#' @param transcript_df Data frame of transcript info. Uses the "cell_id" column to derive cell IDs. 
#' Must provide either this or "cell_ids" argument directly.
#' @export
initPolys <- function(cell_ids = NULL, transcript_df = NULL) {
  if (is.null(cell_ids)) {
    cell_ids <- unique(transcript_df$cell_ids)
  }
  setNames(vector("list", length(cell_ids)), cell_ids)
}

#' Get polygons for specified cell IDs
#' 
#' @param polys Polygons list to be added to
#' @param cell_ids Vector of cell IDs. If NULL, will calc polys for all cells in transcript_df
#' @param transcript_df Data frame of transcript info. Needs these columns: "cell_id", "x", "y".
#' @param overwrite Logical, for whether to overwrite polygons that have already been calculated
#' @param type Either "ahull" or "chull". The latter is faster but limited to convex shapes.
#' @param alpha Passed to alphahull::ahull(). Smaller makes for more nuanced polygons.
#' @return A list of polygons with the entries for the cell IDs filled in
#' @export
cellPolys <- function(polys, cell_ids = NULL, transcript_df, overwrite = FALSE, type = "chull", alpha = 0.1) {
  
  if (is.null(cell_ids)) {
    cell_ids <- unique(transcript_df$cell_id)
  }
  for (cell in cell_ids) {
    if (is.null(polys[[cell]]) || overwrite) {
      use <- transcript_df$cell_id == cell
      polys[[cell]] <- getPolyFromTranscripts(x = transcript_df[use, "x"], 
                                              y = transcript_df[use, "y"], 
                                              type = type, alpha = alpha) 
      
    }
  }
  return(polys)
}


