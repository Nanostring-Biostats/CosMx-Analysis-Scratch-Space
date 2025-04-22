
#' Make a matrix of scores from metagenes object
#' @param metagenes metagenes object
#' @param to_model variable in metagenes object used to make matrix of scores, i.e., "yhat" or "ypost"
#' 
#' @export  
make_scores_matrix <- function(metagenes, to_model){
  scores <- as.matrix(do.call(cbind, lapply(metagenes, "[[", to_model)))
  colnames(scores) <- names(metagenes)
  return(scores)
}