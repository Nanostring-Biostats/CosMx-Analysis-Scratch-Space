
#' Identify genes useful for supervised classification of closely-related cell types 
#' @param mat Expression matrix from just the cell type in question, cells in rows and genes in columns.
#' @param ref Reference matrix for the cell types in question
#' @param xy Spatial coordinates of cells  
#' @return A vector of gene name for use in supervised cell subtyping
getSubtypingGenes <- function(mat, ref, xy) {
  
} 






#' Identify genes useful for unsupervised clustering of closely-related cell types 
#' Intended for subclustering of a broad cell type, e.g. of immune cells into T-cells, B-cells, etc... 
#' @param mat Expression matrix from just the cell type in question, cells in rows and genes in columns.
#' @param xy Spatial coordinates of cells 
#' @return A vector of gene named for use in unsupervised subclustering
getSubclusteringGenes <- function(mat, xy) {
  
} 


#' 