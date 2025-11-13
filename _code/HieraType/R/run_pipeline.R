

#' Run a hierarchical cell typing pipeline 
#' 
#' @param pipeline object created using `make_pipeline()`
#' @param counts_matrix a cells x genes expression matrix.  
#'                  If modeling scaled pearson residuals (pearson.response = TRUE  by default)
#'                  ,this should be a raw counts matrix.
#' @param adjacency_matrix an optional cells x cells matrix matrix of weights denoting similarity between pairs of cells.
#'                         For example, one could use the graph of smoothed nearest neighbors distances used to optimize UMAP embeddings.
#' @param initial_prior_weights an optional vector of weights (typically in range of 0-1), 
#' denoting confidence that the cells belong to any of the classes in provided markerslist.
#' This argument can be used to specify weights for the **first** celltyping step.
#' @param celltype_call_threshold Posterior probability threshold used to determine level of calling for `celltype_thresh` annotation in the returned `post_probs` data.table.
#' See `combine_postprob_tables()` function; these tables can be remade quickly at different thresholds without the need to rerun a pipeline.
#' @param return_all_columns_postprobs Whether to return individual posterior probabilities for every celltype in the `post_probs` data.table. 
#' See `combine_postprob_tables()` function; these columns can be returned later if needed without needing to rerun a pipeline.
#' @param ... Other arguments to be passed to `fit_metagene_scores` and/or `cluster_metagenes` functions.
#' 
#' @export
run_pipeline <- function(pipeline, counts_matrix
                         ,adjacency_matrix = NULL
                         ,initial_prior_weights = NULL
                         ,celltype_call_threshold = 0.5
                         ,return_all_columns_postprobs = FALSE
                         , ...){
  stopifnot("`pipeline` must have class 'pipeline', typically creaed with `make_pipeline()` function. " = inherits(pipeline, "pipeline"))
 
  dots <- list(...) 


   
  ### markerslists which dont inherit from another markerslist are the starting points
  parent_lists <- setdiff(names(pipeline$markerslists), names(pipeline$priors))
  metagene_scores <- models <- vector(mode = 'list',length=length(pipeline$markerslists)) 
  names(metagene_scores) <- names(models) <- names(pipeline$markerslists)
  child_categories <- c()
  for(parnt in parent_lists){
     metagene_scores[[parnt]]  <- 
         do.call(fit_metagene_scores
                 ,c(list(markerslist = pipeline$markerslist[[parnt]]
                         ,counts_matrix = counts_matrix
                         ,adjacency_matrix = adjacency_matrix
                         ,prior_level_weights = initial_prior_weights
                         )
                    ,dots[names(dots) %in% names(formals(fit_metagene_scores))]
                    )
                 )
     models[[parnt]] <- 
       do.call(cluster_metagenes
               ,c(list(
                  metagenes = metagene_scores[[parnt]]
                  ,prior_prob_level = initial_prior_weights
                  )
                  ,dots[names(dots) %in% names(formals(cluster_metagenes))]
                  )
               ) 
       
     child_categories <- c(child_categories, names(pipeline$priors)[which(pipeline$priors == parnt)])
  }
 
  #remaining_catg <- setdiff(names(pipeline$markerslists), names(models)) 
  while(length(child_categories) > 0){
    child_categories_new <- c()
    for(chld in child_categories){
      prior_wts <- models[[pipeline$priors[[chld]]]]$post_probs[[pipeline$priors_category[[chld]]]]
       metagene_scores[[chld]]  <-
         do.call(fit_metagene_scores
                 ,c(list(markerslist = pipeline$markerslist[[chld]]
                         ,counts_matrix = counts_matrix
                         ,adjacency_matrix = adjacency_matrix
                         ,prior_level_weights = prior_wts
                         )
                    ,dots[names(dots) %in% names(formals(fit_metagene_scores))]
                    )
                 )
       
       models[[chld]] <- 
         do.call(cluster_metagenes
                 ,c(list(
                    metagenes = metagene_scores[[chld]]
                    ,prior_prob_level = prior_wts
                    )
                    ,dots[names(dots) %in% names(formals(cluster_metagenes))]
                    )
                 ) 
           
       child_categories_new <- c(child_categories_new, names(pipeline$priors)[which(pipeline$priors == chld)])
    }
    child_categories <- child_categories_new
  }

   
  ### Make output posterior probability tables
  post_probsl <- 
  combine_postprob_tables(pipeline = pipeline
                          ,models = models
                          ,return_all_columns = return_all_columns_postprobs
                          ,celltype_call_threshold = celltype_call_threshold
                          )
  return(list(post_probs = post_probsl
              ,models = models
              ,metagene_scores = metagene_scores))
}


