

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
#' @param level_args Named list of level-specific arguments. Names should match pipeline level names.
#'        Each element should be a list of arguments to override for that level.
#'        These override the global defaults provided via `...`.
#'        Example: `level_args = list(l1 = list(fit_single_positive_only = TRUE))`
#' @param ... Other arguments to be passed to `fit_metagene_scores` and/or `cluster_metagenes` functions.
#'        These serve as global defaults for all pipeline levels. Use `level_args` to override for specific levels.
#'
#' @return A list containing:
#' \itemize{
#'   \item \code{post_probs} - Combined posterior probability tables from all pipeline stages
#'   \item \code{models} - List of fitted clustering models for each pipeline stage
#'   \item \code{metagene_scores} - List of metagene score fits for each pipeline stage
#' }
#'
#' @examples
#' \dontrun{
#' # Basic usage with global defaults
#' result <- run_pipeline(pipeline, counts, adjacency,
#'                        k_components = 6)
#'
#' # Level-specific arguments
#' result <- run_pipeline(pipeline, counts, adjacency,
#'                        k_components = 6,  # default for all levels
#'                        level_args = list(
#'                          l1 = list(fit_single_positive_only = TRUE, k_components = 8),
#'                          tcellmajor = list(k_components = 4)
#'                        ))
#' }
#'
#' @export
run_pipeline <- function(pipeline, counts_matrix
                         ,adjacency_matrix = NULL
                         ,initial_prior_weights = NULL
                         ,celltype_call_threshold = 0.5
                         ,return_all_columns_postprobs = FALSE
                         ,level_args = list()
                         , ...){
  stopifnot("`pipeline` must have class 'pipeline', typically created with `make_pipeline()` function. " = inherits(pipeline, "pipeline"))

  dots <- list(...)

  # Warn on unrecognized ... arguments
  recognized_args <- unique(c(
    names(formals(fit_metagene_scores)),
    names(formals(cluster_metagenes))
  ))
  dropped <- setdiff(names(dots), recognized_args)
  if (length(dropped) > 0) {
    warning("Unrecognized arguments in `...` will be ignored: ",
            paste(dropped, collapse = ", "),
            "\nRecognized arguments are forwarded to fit_metagene_scores() and cluster_metagenes().",
            immediate. = TRUE)
  }

  # Validate level_args names
  if (length(level_args) > 0) {
    invalid_levels <- setdiff(names(level_args), names(pipeline$markerslists))
    if (length(invalid_levels) > 0) {
      warning(paste0("level_args contains unknown pipeline levels: ",
                     paste(invalid_levels, collapse = ", "),
                     ". These will be ignored."),
              immediate. = TRUE)
    }
  }

  # Validate and align initial_prior_weights with counts_matrix
  if (!is.null(initial_prior_weights)) {
    if (!is.null(names(initial_prior_weights)) && !is.null(rownames(counts_matrix))) {
      if (!all(rownames(counts_matrix) %in% names(initial_prior_weights))) {
        stop("initial_prior_weights is missing entries for some cells in counts_matrix")
      }
      initial_prior_weights <- initial_prior_weights[rownames(counts_matrix)]
    } else if (length(initial_prior_weights) != nrow(counts_matrix)) {
      stop("length of initial_prior_weights (", length(initial_prior_weights),
           ") does not match nrow(counts_matrix) (", nrow(counts_matrix), ")")
    }
  }

  # Helper function to get merged arguments for a specific level
  # Level-specific args override global defaults
  get_level_args <- function(level_name, target_fun) {
    valid_args <- names(formals(target_fun))
    args <- dots[names(dots) %in% valid_args]
    if (level_name %in% names(level_args)) {
      level_specific <- level_args[[level_name]]
      level_specific <- level_specific[names(level_specific) %in% valid_args]
      args[names(level_specific)] <- level_specific
    }
    return(args)
  }

  ### markerslists which dont inherit from another markerslist are the starting points
  parent_lists <- setdiff(names(pipeline$markerslists), names(pipeline$priors))
  metagene_scores <- models <- vector(mode = 'list',length=length(pipeline$markerslists))
  names(metagene_scores) <- names(models) <- names(pipeline$markerslists)
  child_categories <- c()
  for(parnt in parent_lists){
     metagene_scores[[parnt]]  <-
         do.call(fit_metagene_scores
                 ,c(list(markerslist = pipeline$markerslists[[parnt]]
                         ,counts_matrix = counts_matrix
                         ,adjacency_matrix = adjacency_matrix
                         ,prior_level_weights = initial_prior_weights
                         )
                    ,get_level_args(parnt, fit_metagene_scores)
                    )
                 )
     models[[parnt]] <-
       do.call(cluster_metagenes
               ,c(list(
                  metagenes = metagene_scores[[parnt]]
                  ,prior_prob_level = initial_prior_weights
                  )
                  ,get_level_args(parnt, cluster_metagenes)
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
                 ,c(list(markerslist = pipeline$markerslists[[chld]]
                         ,counts_matrix = counts_matrix
                         ,adjacency_matrix = adjacency_matrix
                         ,prior_level_weights = prior_wts
                         )
                    ,get_level_args(chld, fit_metagene_scores)
                    )
                 )

       models[[chld]] <-
         do.call(cluster_metagenes
                 ,c(list(
                    metagenes = metagene_scores[[chld]]
                    ,prior_prob_level = prior_wts
                    )
                    ,get_level_args(chld, cluster_metagenes)
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


