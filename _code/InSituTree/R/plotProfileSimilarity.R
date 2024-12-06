#' Plots a heatmap of cosine similarity between cell types
#' in a reference profile
#'
#' @param full_profiles Cell profile matrix, genes x cell types
#' 
#' @export
#'
#' @examples
#' library(InSituType)
#' data("ioprofiles")
#' plotProfileSimilarity(ioprofiles)
#'
plotProfileSimilarity <- function(full_profiles){
  if(!is.matrix(full_profiles)){
    stop("Function expects a matrix.")
  }
  
  ## Cosine similarity across all pairs of columns
  
  # Initialize matrix
  n <- ncol(full_profiles)
  cosine_sim_matrix <- matrix(NA, n, n)
  
  # Calculate cosine similarity
  for (i in 1:n) {
    for (j in 1:n) {
      cosine_sim_matrix[i, j] <- cosine_similarity(full_profiles[, i], full_profiles[, j])
    }
  }
  
  # Add names
  colnames(cosine_sim_matrix) = rownames(cosine_sim_matrix) = colnames(full_profiles)
  
  # Plot it
  heatmap_res <- heatmap(cosine_sim_matrix)
  return(heatmap_res)
}

#' Calculate cosine similarity between two vectors
#' 
#' @param x vector 1
#' @param y vector 2
#'
#' @export
#' 
cosine_similarity <- function(x, y) {
  sum(x * y) / (sqrt(sum(x^2)) * sqrt(sum(y^2)))
}
