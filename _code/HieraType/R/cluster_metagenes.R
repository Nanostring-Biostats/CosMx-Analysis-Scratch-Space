

#' Cluster metagene scores using a mixture of overfitted gaussian mixture models.
#'
#' @param metagenes metagenes object created by
#' @param scores optional cells x cell types matrix of metagene scores, only needed if `metagenes` not provided.
#' @param to_model either 'yhat' (predicted scores, recommended) or 'ypost' posterior-mean metagene scores.
#' @param discourage_double_positive list in the format of
#' list("mycelltype" = c("celltypeA_shouldnt_be_positive_in_mycelltype", "celltypeB_shouldnt_be_positive_in_mycelltype"))
#' In this example, cells with positive metagene scores for "celltypeA" and/or "celltypeB" will not be used as anchor cells in "mycelltype"
#' This effectively encourages double positive cells to be assigned to "celltypeA" or "celltypeB" in the event they are double-positive with "mycelltype".
#' This can be useful if marker genes for "celltypeA","celltypeB" are weak or low expressors, and/or mycelltype has a high expressing promiscuous marker gene which may bleed into other celltypes.
#' @param prior_prob_level optional vector of weights / posterior probabilities (likely from a previous cluster_metagenes model)
#' that the cells belong to *any* of the classes in 'metagenes'.
#' @param fit_single_positive_only If TRUE, then models will be trained using cells which are *only* positive for that corresponding cell type.
#' Otherwise, cells which are positive in both observed AND fitted values for a particular celltype metagene will be included in model training
#' (even if they have positive scores for other cell types), provided that they pass any `discourage_double_positive` constraints.
#' @param fit_threshold_prior_prob_level If prior_prob_level is provided, this is the minimum threshold considered for a cell to be included in
#' mixture model training for a particular cell type.
#' @param k_components Number of mixture components to be used in fitting mixture model for each cell type.
#' Goal is to have 'enough', to flexibly model the distribution of metagene scores.
#' @param training_cellids  If provided, the metagenes or scores object will be subset to only this set of cells before clustering and scoring.
#' @param seed  seed to be set in local environment before using k-means to initialize starting parameters.  Ensures consistency after re-running.
#' Global seed is reset on function exit.
#'
#' @export
cluster_metagenes <- function(metagenes = NULL
                              ,scores = NULL
                              ,to_model = c("yhat", "ypost")
                              ,discourage_double_positive = list("plasma" = "immune"
                                                                 ,"endothelial" = "immune"
                                                                 ,"epithelial" = "immune"
                                                                 ,"fibroblast" = "immune"
                                                                 ,"smooth_muscle" = "immune"
                                                                 ,"mast" = "tcell")
                              ,prior_prob_level = NULL
                              ,fit_single_positive_only = FALSE
                              ,fit_threshold_prior_prob_level = 0.1
                              ,verbose = 1
                              ,k_components = 6
                              ,training_cellids = NULL
                              ,seed = 123
                              ){

  ### scores to be modeled 
  to_model <- match.arg(to_model) 
  
  orig_seed <- .Random.seed
  on.exit({.Random.seed <<- orig_seed})
  
  
  if(is.null(scores)){
    if(is.null(metagenes)) stop("either 'scores' or 'metagenes' object must be provided.")
    scores <- make_scores_matrix(metagenes, to_model)
  }

  if(!is.null(prior_prob_level)){
    stopifnot(nrow(scores) == length(prior_prob_level))
    if(is.null(names(prior_prob_level))){
      names(prior_prob_level) <- rownames(scores)
    }  else {
      stopifnot("mismatch between metagene cell ids and names provided in prior_prob_level" = 
                  all(rownames(scores) %in% names(prior_prob_level)))
      prior_prob_level <- prior_prob_level[rownames(scores)]
    }
  } 
  
  npos <- rep(1, nrow(scores))   ## holds number of positive metagene scores per cell 
  npos.2 <- rep(1, nrow(scores)) ## holds logical condition whether `prior_prob_level` exceedds the minimum threshold `fit_threshold_prior_prob_level` 
  
  ### identify single-positive anchor cells to fit model to
  if(!is.null(training_cellids)){
    scores <- scores[training_cellids,] 
    npos <- npos[1:length(training_cellids)]
    npos.2 <- npos.2[1:length(training_cellids)]
  } else {
    npos <- apply(scores , 1, function(x) sum(x > 0))
    if(!is.null(prior_prob_level)){
      npos.2 <- as.numeric(prior_prob_level > fit_threshold_prior_prob_level)
    }
  }

  anypos <- Matrix::colSums(scores > 0)
  anchorid <- names(anypos[anypos > 0])
  models <- vector(mode = 'list', length=length(unique(anchorid)))
  names(models) <- intersect(colnames(scores), unique(anchorid))
  nfit <- vector(mode = 'numeric', length=length(models))
  names(nfit) <- names(models)
  yh <- make_scores_matrix(metagenes, "yhat")
  y <- make_scores_matrix(metagenes, "y")
  yp <- make_scores_matrix(metagenes, "ypost")
  conds <- yh > 0 & y > 0 ## both observed and predicted score are positive
  for(ii in intersect(colnames(scores), unique(anchorid))){
    criteria_1 <- (npos.2==1 & npos == 1 & (scores[,ii] > 0))
    scoresfit_anchor <- scores[criteria_1,]
    if(!fit_single_positive_only){
      criteria_2 <- (npos.2==1 & conds[,ii]==1) ## passes prior_prob_level threshold and yhat > 0 and ypost(erior) > 0
      if(!is.null(discourage_double_positive[[ii]])){
        if(any(discourage_double_positive[[ii]] %in% colnames(yp))){
          otherct <- intersect(discourage_double_positive[[ii]],colnames(yp))
          double_pos_pass_yhat <- apply(yp[,otherct,drop=FALSE] > 0, 1, mean) ## check if posterior scores for forbidden cell types are > 0
          double_pos_pass <- double_pos_pass_yhat==0  ## if none of the forbidden celltypes have posterior scores > 0, then the cell passes
          criteria_2 <- criteria_2 & double_pos_pass  
        } else {
          warning(paste0("`discourage_double_positive` argument ignored.  Specified for "
                         , ii, ",\nbut celltypes to discourage ("
                         ,paste0(discourage_double_positive[[ii]], collapse=",")
                         , ") not found."))
        }
      }
      scoresfit_anchor <- scores[criteria_1 | criteria_2,]  ## either single-positive OR [ both yhat AND ypost are positive AND none of the forbidden cell types are positive] 
    }
    nfit[ii] <- nrow(scoresfit_anchor)
    if(nfit[ii] > 0){
      message(paste0("fitting mixture component for ", ii))
      ktry <- k_components; tryagain <- 1
      while(tryagain & ktry > 0){
        tryagain <- 0
        set.seed(seed)
        inits <- kmeans_init(scoresfit_anchor, k = ktry)
        ### check for singletons
        singleton_has_na <- unlist(lapply(inits$initsigma, function(x) sum(is.na(x))))
        if(any(singleton_has_na > 0)){
          inits$initmu <- inits$initmu[-c(which(singleton_has_na > 0)),,drop=FALSE]
          inits$initsigma <- inits$initsigma[-c(which(singleton_has_na > 0))]
        }
        
        models[[ii]] <- 
        em_mvn_mixture_overfit(scoresfit_anchor
                               ,initmu = inits$initmu
                               ,initsigma = inits$initsigma
                               ,prior_prob_level = prior_prob_level[rownames(scoresfit_anchor)]
                               )
        if(any(models[[ii]]$pihat==0)){
          tryagain <- 1
          ktry <- ktry - 1
        }
      }
    }
  }
  
  priork <- prop.table(nfit)
  ll_scores <- vector(mode = 'list', length=length(models))
  names(ll_scores) <- names(models)
  
  for(ii in names(models)){
    ll_scores[[ii]] <- 
    score_mixture_overfit(models[[ii]]
                          ,scores = scores
                          ,denom_only = TRUE ## total loglik across K clusters
                          )
    ll_scores[[ii]] <- ll_scores[[ii]] + log(priork[ii])
  }
  llmat <- do.call(cbind, ll_scores)
  ppmatdenom <- apply(llmat, 1, matrixStats::logSumExp)
  ppmatnum <- exp(llmat)
  ppmat <- Matrix::Diagonal(x=1/exp(ppmatdenom),names=TRUE)%*% ppmatnum  
  
  isna <- apply(ppmat, 1, function(x) sum(is.na(x) | is.infinite(x)))
  if(sum(isna > 0) > 0){
    whichisna <- which(isna > 0)
    maxll <- apply(llmat[whichisna,,drop=FALSE], 1, which.max)
    for(j in 1:length(whichisna)){
      ppmat[whichisna[j],] <- 0
      ppmat[whichisna[j],maxll[j]] <- 1
    }
  } 
  if(!is.null(prior_prob_level)){
    ppmat <- (Matrix::Diagonal(x=prior_prob_level,names=TRUE) %*% ppmat)
  } 
  
  post_probs <- data.table::data.table(as.matrix(ppmat))
  colnames(post_probs) <- colnames(scores)
  post_probs[,best_score:=do.call(pmax,.SD)]
  post_probs[,best_class:=colnames(post_probs)[apply(.SD,1,which.max)],.SDcols=1:(ncol(post_probs)-1)]
#  post_probs[,best_class:=colnames(post_probs)[which.max(.SD)],by=.I,.SDcols=(1:(ncol(post_probs)-1))]
  post_probs[,cell_ID:=rownames(llmat)] 
  data.table::setcolorder(post_probs, c("cell_ID", "best_class", "best_score")) 
  lx_model <- list(models = models, post_probs = post_probs, llmat = llmat) 
  
  ### return results
  return(lx_model) 
}
