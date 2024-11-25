#' Summarize the annotations for runInSituTree() function
#'
#' @param nestList Nested insitutype objects returned from runInSituType()
#' @param annotationLevel Used for recursive looping.  Do not adjust.
#'
#' @return Dataframe of celltype annotation and posterior probabilities
#' @export
#' @import dplyr
#' @import tibble
#'
#'
#'
summarizeInSituTree <- function(nestList, annotationLevel = 1){
  message(paste0("Summarizing group: ", nestList$name))
  df_out <- data.frame(nestList$result$clust, nestList$result$prob)
  colnames(df_out) <- c(paste0("annotLevel_", annotationLevel), paste0("probs_", "annotLevel_", annotationLevel))

  # summarize all subclusterings
  if(length(nestList$subclusterings) >= 1){
    sub_dfs <- lapply(nestList$subclusterings, function(gg){
      message(paste0(unique(gg$result$clust), collapse = " & "))
      summarizeInSituTree(nestList = gg, annotationLevel = annotationLevel + 1)
    })
    sub_df <- Reduce(dplyr::bind_rows, sub_dfs)

    # merge the new celltype results with the higher level annotations.

    # Ensure row names are columns before merging
    df_out <- df_out %>%
      tibble::rownames_to_column(var = "Row.names")
    sub_df <- sub_df %>%
      tibble::rownames_to_column(var = "Row.names")

    # Merge data frames using dplyr's full_join to keep all rows
    df_out <- full_join(df_out, sub_df, by = "Row.names")

    # Remove the row names column after merging
    df_out <- df_out %>%
      tibble::column_to_rownames(var = "Row.names")

    # copy over the higher level annotations for cells that only have higher level annotations.
    for(i in 2:(ncol(df_out)/2)){
      if(any(is.na(df_out[,i*2-1]))){
        df_out[is.na(df_out[,i*2-1]),i*2-1] <- df_out[is.na(df_out[,i*2-1]),i*2-3]
        df_out[is.na(df_out[,i*2]),i*2] <- df_out[is.na(df_out[,i*2]),i*2-2]
      }
    }
  }



message(paste0("Completed summarizing group: ", nestList$name))

return(df_out)

}
