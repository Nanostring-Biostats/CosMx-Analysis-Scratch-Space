
#' Utility function for making / remaking or restyling 
#' combined posterior probability tables from clustering models (possibly at different thresholds).
#' 
#' @param pipeline a pipeline object created using `make_pipeline()`
#' @param models list of models, created using `run_pipeline`
#' @param celltype_call_threshold Used to determing the celltype_thresh column in returned table.  
#' This column shows the most granular cell type with posterior probability > `celltype_call_threshold` for each cell.
#' @param return_all_columns if FALSE (default), columns with posterior probability scores for
#' individual cell type classes are dropped from the returned table.
#' @return A named list of data.tables (one per parent pipeline), each containing cell_ID,
#'         celltype_thresh, celltype_granular, best_score_thresh, and best_score_granular columns.
#' @export
combine_postprob_tables <- function(pipeline
                                    ,models
                                    ,celltype_call_threshold = 0.5
                                    ,return_all_columns = FALSE
                                    ){
  stopifnot("`pipeline` must have class 'pipeline', typically creaed with `make_pipeline()` function. " = inherits(pipeline, "pipeline"))
  parent_lists <- setdiff(names(pipeline$markerslists), names(pipeline$priors))
  
  ### Make output posterior probability tables
  post_probsl <- vector(mode='list',length=length(parent_lists))
  names(post_probsl) <- parent_lists
  for(parnt in parent_lists){
    post_probsl[[parnt]] <- copy(models[[parnt]]$post_probs) 
    post_probsl[[parnt]][,celltype_thresh:=parnt]
    post_probsl[[parnt]][,celltype_granular:=best_class]
    post_probsl[[parnt]][best_score > celltype_call_threshold,celltype_thresh:=best_class]
    post_probsl[[parnt]][,best_score_granular:=best_score]
    post_probsl[[parnt]][,best_score_thresh:=best_score]
    #post_probsl[[parnt]][best_score_thresh < celltype_call_threshold,best_score_thresh:=NA_real_]
    data.table::setnames(post_probsl[[parnt]], c("best_score", "best_class"), paste0(c("best_score_", "best_class_"), parnt))
    child_categories <- c(names(pipeline$priors)[which(pipeline$priors == parnt)])
    while(length(child_categories) > 0){
      child_categories_new <- c()
      for(chld in child_categories){
        post_probsl[[parnt]] <- merge(post_probsl[[parnt]], models[[chld]]$post_probs, by = "cell_ID", sort = FALSE)
        #post_probsl[[parnt]][best_score > celltype_call_threshold,celltype_thresh:=best_class]
        whch_threshold <- post_probsl[[parnt]][best_score > celltype_call_threshold,which=TRUE]
        parnt_cat <- pipeline$priors_category[[chld]]
        parnt_lst <- pipeline$priors[[chld]]
        whch <- which(post_probsl[[parnt]][[paste0("best_class_",parnt_lst)]] == parnt_cat)
        
        go_up_tree_n <- 1
        while(go_up_tree_n > 0){
          if(!is.null(pipeline$priors[[parnt_lst]])){
            parnt_cat <- pipeline$priors_category[[parnt_lst]]
            parnt_lst <- pipeline$priors[[parnt_lst]]
            whch <- intersect(whch, which(post_probsl[[parnt]][[paste0("best_class_",parnt_lst)]] == parnt_cat))
            whch_threshold <- intersect(whch_threshold, post_probsl[[parnt]][best_score > celltype_call_threshold,which=TRUE])
          } else {
            go_up_tree_n <- 0
          }
        }
        post_probsl[[parnt]][whch,celltype_granular:=best_class]
        post_probsl[[parnt]][whch,best_score_granular:=best_score]
        post_probsl[[parnt]][intersect(whch_threshold, whch),celltype_thresh:=best_class]
        post_probsl[[parnt]][intersect(whch_threshold, whch),best_score_thresh:=best_score]
        data.table::setnames(post_probsl[[parnt]], c("best_score", "best_class"), paste0(c("best_score_", "best_class_"), chld))
        child_categories_new <- c(child_categories_new, names(pipeline$priors)[which(pipeline$priors == chld)])
        if(!return_all_columns){
          keepcols <- grep("^celltype_|^best_score_|^best_class_|cell_ID", names(post_probsl[[parnt]]), value = TRUE)
          post_probsl[[parnt]] <- post_probsl[[parnt]][,keepcols,with=FALSE]
        }
      } 
      child_categories <-  child_categories_new
    }
    if(!return_all_columns){
      keepcols <- grep("^celltype_|^best_score_thresh|^best_score_granular|cell_ID", names(post_probsl[[parnt]]), value = TRUE)
      post_probsl[[parnt]] <- post_probsl[[parnt]][,keepcols,with=FALSE]
    }
  } 
  return(post_probsl)
}

