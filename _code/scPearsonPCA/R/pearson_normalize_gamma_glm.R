
#' Normalize continuous protein expression data using pearson residuals
#' 
#' @description
#' Fits a gamma generalized linear model for each protein / row in the expression matrix `x`, 
#' and returns the pearson residuals.  This new proteins x cells matrix can be used as a normalized dataset for 
#' downstream analyses.
#'
#' @param x A proteins x cells expression matrix 
#' @param lower_q_thresh  lower quantile at which to threshold raw expression for a given protein (default = 0.01) before fitting gamma glm model.  
#' Note that for a gamma family glm used here, the protein expression must be positive.
#' @param upper_q_thresh  upper quantile at which to threshold expression for a given protein (none by default) before fitting gamma glm model.  
#' To avoid outliers skewing normalization, could threshold extreme values of raw protein expression
#'
#'
#' @export  
#' 
pearson_normalize_gamma_glm <- function(x
                                        ,lower_q_thresh = 0.01
                                        ,upper_q_thresh = NULL
){
  
  if(is.null(lower_q_thresh)) lower_q_thresh <- 0 
  if(is.null(upper_q_thresh)) upper_q_thresh <- 1
  stopifnot(upper_q_thresh <= 1)
  stopifnot(lower_q_thresh >= 0)
  
  ## must prevent <=0 for Gamma model
  lower_q <- apply(x, 1, quantile, lower_q_thresh)
  lower_q[lower_q < 0] <- min(lower_q[lower_q > 0])
  
  if(upper_q_thresh < 1){
    upper_q <- apply(x, 1, quantile, upper_q_thresh)
  } 
  
  pgamm_res <- 
    lapply(1:nrow(x), function(ii){
      message(ii)
      ymodel <- x[ii,]
      ymodel <- pmax(ymodel, lower_q[ii]) ## must prevent <=0
      if(upper_q_thresh < 1){  ## optionally trim to upper quantile
        ymodel <- pmin(ymodel, upper_q[ii])
      }
      gm <- glm(ymodel ~ offset(log(Matrix::colSums(x)))
                ,family = Gamma(link = "log"))
      return(residuals(gm,type="pearson"))
    })
  
  pgamm_res <- do.call(rbind, pgamm_res)
  gc()
  rownames(pgamm_res) <- rownames(x)
  colnames(pgamm_res) <- colnames(x)
  return(pgamm_res) 
}

