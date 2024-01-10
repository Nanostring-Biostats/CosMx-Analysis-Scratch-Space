
#' Update cell typing results with spatial context or other alternative data
#' 
#' Takes cell typing results, then updates it based on alternative data types, 
#' e.g. spatial context, morphology, or proteomics. Existing cell typing results are 
#' put into Insitutype's likelihood framework, which then can use alternative data
#' as a prior to be updated by the expression data to get a new posterior probability 
#' of cell type.
#' Performs this operation by
#' 1. deriving cell type profiles using InSituType:::Estep(), 
#' 2. assigning cells to "cohorts" (clusters) derived from their alternative data
#' 3. Inputing the output of steps (1) and (2) into InSituType::insitutype() to 
#'  re-calculate cell type. 
#' @param
#' 
spatialUpdate <- function() {
  
}