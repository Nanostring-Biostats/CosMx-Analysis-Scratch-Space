#' Make a stained glass plot
#' 
#' Uses drawPolys to make a plot resembling stained glass
#' @param polys List of polygons
#' @param transcript_df Data frame of transcript info. Needs these columns: "cell_id", "x", "y".
#' @param celltype Vector of cell types, with names aligned to the names of polys
#' @param transcript_cex Point size of transcripts. Default works well for a 1/4 FOV view
#' @param cellcols Vector of colors, named by celltype
stainedGlassPlot <- function(polys, transcript_df, celltype, 
                             transcript_cex = 0.3, cellcols = notredamecols) {
  
  # set colors:
  ndcellcols = notredamecols[1 + (1:length(unique(celltype))) %% length(notredamecols)]
  names(ndcellcols) = names(table(celltype))[order(table(celltype), decreasing = T)]
  
  bgpar <- par()$bg
  par(bg = "black")
  par(mar = c(0,0,0,0))
  plot(transcript_df[, c("x", "y")],
       pch = 16, col = scales::alpha("white", 0.1), cex = transcript_cex, asp = 1)
  drawPolys(polys, col = scales::alpha(ndcellcols[celltype[names(polys)]], 0.6), 
            border = "black", 
            add = TRUE, asp = 1)
  par(bg = bgpar)
}

