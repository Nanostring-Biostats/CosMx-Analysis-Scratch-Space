

#' check that genes are available and remove them if missing.
#'  
validate_markerslist <- function(markerslist, available_genes){
 
  ## check for missing index markers.
  index_markers <- unlist(lapply(markerslist, "[[", "index_marker"))
  missing_index_markers <- setdiff(index_markers, available_genes)
  if(length(missing_index_markers) > 0){
    stop(paste0("The specified index markers: ", paste0(missing_index_markers, collapse=","), " were not found."))
  }
  
  ## check for missing predictors.  Remove any missing and print a message.
  all_missing_predictors <- c()
  for(ct in names(markerslist)){
      predct <- markerslist[[ct]][["predictors"]] 
      not_missing_predictors <- which(predct %in% available_genes)
      if(length(not_missing_predictors) < length(predct)){
        markerslist[[ct]][["predictors"]] <- predct[not_missing_predictors]
        markerslist[[ct]][["signs"]] <- markerslist[[ct]][["signs"]][not_missing_predictors]
        all_missing_predictors <- unique(c(all_missing_predictors, setdiff(predct, predct[not_missing_predictors])))
      }
  }
  
  if(length(all_missing_predictors) > 0){
      message(paste0("Removed ", length(all_missing_predictors), " predictor genes not found in the counts matrix: "
              ,paste0(all_missing_predictors, collapse=", ")))
  }

  positive_predictor_lengths <- 
  unlist(lapply(markerslist, function(xx){
    sum(xx$signs==1)
  }))
  if(any(positive_predictor_lengths < 5)){
    for(ii in names(positive_predictor_lengths[positive_predictor_lengths < 5])){
      warning(paste0(ii, " has only ", positive_predictor_lengths[ii], " positive predictor genes.  Consider adding more predictor genes for this celltype.")) 
    } 
  }
  return(markerslist) 
}






