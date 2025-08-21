

#' check that genes are available and remove them if missing.
#'  
validate_markerslist <- function(markerslist, available_genes){

  ## check for missing index markers.
  index_markers <- lapply(markerslist, "[[", "index_marker")
  for(idx in 1:length(index_markers)){
    missing_index_markers <- setdiff(index_markers[[idx]], available_genes)
    if(length(missing_index_markers) > 0){
      if(length(missing_index_markers)==length(index_markers[[idx]])){
        stop(paste0(names(index_markers)[idx], ": All specified index markers: ", paste0(missing_index_markers, collapse=","), " were not found."))
      } else {
        markerslist[[idx]][["index_marker"]] <- setdiff(index_markers[[idx]], missing_index_markers)
       # markerslist[[idx]][["predictors"]] <- markerslist[[idx]][["predictors"]][c(setdiff(index_markers[[idx]], missing_index_markers))]
#        markerslist[[idx]][["signs"]] <- markerslist[[idx]][["signs"]][c(setdiff(index_markers[[idx]], missing_index_markers))]
        msg <- paste0(names(index_markers)[idx], ":  the specified index markers: ", paste0(missing_index_markers, collapse=","), " were not found.")
        msg <- paste0(msg, "\n", length(markerslist[[idx]][["index_marker"]]), " remaining:", paste0(markerslist[[idx]][["index_marker"]], collapse=","))
        warning(msg)
      }
      
    }
  } 
  ## check for missing predictors.  Remove any missing and print a message.
  all_missing_predictors <- c()
  for(ct in names(markerslist)){
    predct <- markerslist[[ct]][["predictors"]]
    not_missing_predictors <- which(predct %in% available_genes)
#   or(idx in 1:length(markerslist[[ct]][["predictors"]])){
    if(length(not_missing_predictors) < length(predct)){
      markerslist[[ct]][["predictors"]] <- predct[not_missing_predictors]
     # markerslist[[ct]][["signs"]][[idx]] <- markerslist[[ct]][["signs"]][[idx]][not_missing_predictors]
      all_missing_predictors <- unique(c(all_missing_predictors, setdiff(predct, predct[not_missing_predictors])))
    }
  }
  
  if(length(all_missing_predictors) > 0){
      message(paste0("Removed ", length(all_missing_predictors), " predictor genes not found in the counts matrix: "
              ,paste0(all_missing_predictors, collapse=", ")))
  }

  positive_predictor_lengths <- 
    unlist(lapply(markerslist, function(xx){ length(xx$predictors)}))
  
  if(any(positive_predictor_lengths < 5)){
    for(ii in names(positive_predictor_lengths[positive_predictor_lengths < 5])){
      warning(paste0(ii, " has only ", positive_predictor_lengths[ii], " positive predictor genes.  Consider adding more predictor genes for this celltype.")) 
    } 
  }
  return(markerslist) 
}






