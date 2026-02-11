
#' Create a 'pipeline' object that can be used to chain together hierarchical cell typing steps. 
#' 
#' @param markerslists named list of markerslists objects (created using `make_markerslist()` function), 
#' for each celltyping stage.
#' @param priors if a celltyping stage 'A' depends on another celltyping stage 'B', then specify it as a prior,
#'  i.e., `priors = list("A" = "B")`
#'
#' @param priors_category if a celltyping stage 'A' depends on another celltyping stage 'B', then specify the parent category in 'B' that 'A' is inheriting from ,
#'  i.e., `priors_category = list("A" = "B_type")`
#'
#' @examples
#' 
#' pipeline_tcell <-
#' make_pipeline(markerslists = list("tmajor" = HieraType::markerslist_tcellmajor
#'                                   ,"t4minor" = HieraType::markerslist_cd4tminor
#'                                   ,"t8minor" = HieraType::markerslist_cd8tminor)
#'               ,priors = list("t4minor" = "tmajor"
#'                              ,"t8minor" = "tmajor")
#'               ,priors_category = list("t4minor" = "cd4t"
#'                                       ,"t8minor" = "cd8t")
#' )
#'
#' @return A 'pipeline' object (list with class "pipeline") containing markerslists, priors, and priors_category.
#' @export
#'
#' 
make_pipeline <- function(markerslists
                          ,priors = NULL
                          ,priors_category = NULL
                          ){
  stopifnot(all(!is.null(names(markerslists))))
  if(!is.null(priors)){
    stopifnot("`priors` must be NULL, or a named list, with names corresponding to elements in the `markerslists`" = 
              all(!is.null(names(priors)))
       )
    stopifnot("Names of celltypes in `priors` do not all correspond to names in `markerslists`" =  
                all(names(priors) %in% names(markerslists)))
    
    stopifnot("`priors_category` should be a named list with the same length and same names as `priors`" =  
                all(names(priors_category) %in% names(priors)))
    stopifnot("`priors_category` should be a named list with the same length and same names as `priors`" =  
                all(names(priors) %in% names(priors_category)))
    
    for(celltype in names(priors_category)){
      if(length(priors_category[[celltype]])!=1){
        msg <- paste0("`priors_category` for celltype=",celltype," does not have length 1.\n"
                      ,"each celltype specified with a 'prior' should should inherit from 1 and only "
                      ,"1 category in the 'prior' markerslist.\n"
                      ,"See help(make_pipeline) for examples.")
        stop(msg)
      }
      if(!priors_category[[celltype]] %in% names(markerslists[[priors[[celltype]]]])){
        msg <- paste0("`priors_category` '",priors_category[[celltype]]
                      ,"' not found as a celltype in the prior markerslist '", priors[[celltype]], "'.")
        stop(msg)
      }
    } 
     
    if(length(priors) == length(markerslists)){
      if(all(unlist(lapply(priors, "length")) > 0)){
        stop(paste0("All celltypes in 'markerslists' are specified with a 'prior'.\n"
                    ,"At least one celltype should be omitted from priors, in order to be used as "
                    ,"a starting point for the pipeline."
                    ))
      }
    }
  }
 
   
  for(celltype in names(markerslists)){
    stopifnot(inherits(markerslists[[celltype]], "markerslist"))
  }
  
  pipeline <- list(markerslists = markerslists, priors = priors, priors_category = priors_category)
  class(pipeline) <- append(class(pipeline), "pipeline")
  return(pipeline)
}
