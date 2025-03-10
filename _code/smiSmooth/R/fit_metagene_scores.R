




#' Construct celltype metagene scores for clustering using `make_markerslist` object via non-negative least squares.
#' 
#' @param markerslist a markerslist created using the `make_markerslist` function
#' @param counts_matrix a cells x genes expression matrix.  
#'                      If modeling scaled pearson residuals (pearson.response = TRUE  by default)
#'                      ,this should be a raw counts matrix.
#' @param adjacency_matrix an optional cells x cells matrix matrix of weights denoting similarity between pairs of cells.
#'                         For example, one could use the graph of smoothed nearest neighbors distances used to optimize UMAP embeddings.
#' @param prior_level_weights an optional vector of weights (typically in range of 0-1), denoting confidence that the cells belong to any of the classes in provided markerslist.
#'                            For example, if I am constructing metagene scores for different types of 'tcells', 
#'                            the weights may be posterior probabilities that the cells are 'tcells', which are generated in a previous modeling stage.
#' @param modeltype Right now, just 'nnls' (non-negative least squares) is supported.
#' @param lag_response Should the (weighted) average of the index gene across similar cells be included as a predictor in the regression model? default is TRUE
#' @param lag_predictors Should the (weighted) averages of predictor genes across similar cells be included as predictors in the regression model? default is TRUE
#' @param pearson.response Should the response variable be the pearson residuals of the index variable?
#' @param pearson.scale.max If `pearson.response = TRUE`, have large pearson residuals > `pearson.scale.max` SDs above the mean truncated to this amount.
#' @param gene_wise_frequency If passing in a subset of all genes in the `counts_matrix`, this can be precomputed gene-wise frequencies used to construct pearson residuals.  
#' See `gene_frequency()` for examples in generating this.
#' @param totalcounts If passing in a subset of all genes in the `counts_matrix`, this can be a precomputed vector of 'totalcounts' for each cell, which would be used in calculation of pearson residuals and totalcount normalization. 
#' @param normalize_predictors_by_totalcounts default is TRUE.  
#' @param row_standardize_adjacency_matrix  default is TRUE; row-standardize the adjacency matrix , such that the sum of neighbor weights across rows adds up to 1.
#' @param obs Optional data.frame or data.table of metadata with a column for `batch_variable` and the `cell id`.  If provided and `pearson.response=TRUE` this is used to construct the pearson residuals
#' @param batch_variable Optional column name of 'batch_variable' in `obs` dataframe. If provided and `pearson.response=TRUE` this is used to construct the pearson residuals
#' @param cellid_colname column name of 'cell ids' in the `obs` dataframe.  
#' @return a 'markerslist' list object which can be passed to `fit_metagene_scores`
#'
#' @export
fit_metagene_scores <- function(
                         markerslist
                         ,counts_matrix
                         ,adjacency_matrix = NULL
                         ,prior_level_weights = NULL
                         ,modeltype = c("nnls")
                         ,lag_response = TRUE
                         ,lag_predictors = TRUE
                         ,pearson.response = TRUE
                         ,pearson.scale.max = 10 
                         ,pearson.family = c("poisson", "nonzero_bernoulli")
                         ,totalcounts = NULL
                         ,gene_wise_frequency = NULL
                         ,normalize_predictors_by_totalcounts = TRUE
                         ,row_standardize_adjacency_matrix = TRUE
                         ,obs = NULL
                         ,batch_variable=NULL
                         ,cellid_colname="cell_ID"
                         ,verbose = TRUE
  ){
 
  # check counts in markerslist 
  pearson.family <- match.arg(pearson.family , c("poisson", "nonzero_bernoulli")) 
  
  # prior_level_weights 
  allgenes <- 
  unique(c(unlist(lapply(markerslist, "[[", "predictors"))
          ,unlist(lapply(markerslist, "[[", "index_marker"))))
 
  # check_valid_markerslist / check that all genes present 
  if(!all(allgenes %in% colnames(counts_matrix))){
     if(!all(allgenes %in% rownames(counts_matrix))){
       available_genes <- colnames(counts_matrix) 
       markerslist <- validate_markerslist(markerslist, available_genes)  
       allgenes <- intersect(allgenes, available_genes)
     } else {
       counts_matrix <- Matrix::t(counts_matrix)
       available_genes <- colnames(counts_matrix) 
       markerslist <- validate_markerslist(markerslist, available_genes)  
       allgenes <- intersect(allgenes, available_genes)
     }
  }  else {
    available_genes <- colnames(counts_matrix) 
    markerslist <- validate_markerslist(markerslist, available_genes)  
    allgenes <- intersect(allgenes, available_genes)
  }
 
  if(is.null(totalcounts)){
    totalcounts <- Matrix::rowSums(counts_matrix) 
  }
  stopifnot(length(totalcounts) == nrow(counts_matrix))
  if(!is.null(names(totalcounts)) && !is.null(rownames(counts_matrix))){
    stopifnot("ids of 'totalcounts' dont match ids in 'counts_matrix'" = 
                all(names(totalcounts) %in% rownames(counts_matrix))) 
    totalcounts <- totalcounts[rownames(counts_matrix)] 
  }
  
  if(!is.null(prior_level_weights)){
    stopifnot(length(prior_level_weights) == nrow(counts_matrix))
    if(!is.null(names(prior_level_weights)) && !is.null(rownames(counts_matrix))){
      stopifnot("ids of 'prior_level_weights' dont match ids in 'counts_matrix'" = 
                  all(names(prior_level_weights) %in% rownames(counts_matrix))) 
      prior_level_weights <- prior_level_weights[rownames(counts_matrix)] 
    }
  } 
  
  if(is.null(gene_wise_frequency)){
    if(!is.null(batch_variable)){
      md <- data.table::copy(data.table::data.table(obs))
      if(!(cellid_colname) %in% names(md)){
        stopifnot(!is.null(rownames(obs)))
        cellid_colname <- "cell_ID"
        md[[cellid_colname]] <- rownames(obs)
      }
      if(cellid_colname!="cell_ID" & "cell_ID" %in% names(md)) md[["cell_ID"]] <- NULL
      setnames(md, old=cellid_colname, new="cell_ID")
      
      rm(obs); gc()
      md <- md[match(rownames(counts_matrix),cell_ID)]
      
      md[,grpid__:=.GRP,by=c(batch_variable)]
      batch_mat <- Matrix::sparseMatrix(j = 1:nrow(md), i=md[["grpid__"]], x = 1
                                        ,dimnames = list(c(md[,head(.SD, 1),by=grpid__][order(grpid__)][[batch_variable]])
                                                         ,c(md[["cell_ID"]]))
                                                         
                                        )
      ### expression rates by batch (\hat{p})
      grate <- Matrix::Diagonal(x=1/Matrix::rowSums(batch_mat), names = rownames(batch_mat)) %*% batch_mat %*%  counts_matrix 
      grate <- Matrix::Diagonal(x = 1/Matrix::rowSums(grate),names=rownames(batch_mat)) %*% grate 
    } else {
      gene_wise_frequency <- Matrix::colSums(counts_matrix)
      gene_wise_frequency <- gene_wise_frequency / sum(gene_wise_frequency)
    }
  }
  
 
  if(length(allgenes) < ncol(counts_matrix)){
    counts_matrix <- counts_matrix[,allgenes]  
    gc()
  } 
  
  ### response variable
  y <- counts_matrix[,unlist(lapply(markerslist, "[[", "index_marker"))]

  ## check adjacency matrix compatible with counts matrix
  if(!is.null(adjacency_matrix)){
    if(row_standardize_adjacency_matrix){
      adjacency_matrix <- Matrix::Diagonal(x = pmax(0, 1/Matrix::rowSums(adjacency_matrix))
                                           ,names = rownames(adjacency_matrix)
                                           ) %*% adjacency_matrix
      
    }
     
  }
  
  metagene_fits <- vector(mode='list', length=length(markerslist))
  names(metagene_fits) <- names(markerslist)
  for(jj in 1:length(markerslist)){
    ctclass <- names(markerslist)[jj]
    index_marker <- markerslist[[jj]]$index_marker
    predictors <- markerslist[[jj]]$predictors
    message(ctclass)
    yresp <- counts_matrix[,index_marker]
    if(pearson.response){
      ### summary statistics for normalization and 
      ### (optional) pearson standardization of response variable used for fitting metagene
      if(is.null(batch_variable)){
        muhat <- gene_wise_frequency[index_marker] * totalcounts
      } else {
        muhat <- (Matrix::Diagonal(x = totalcounts, names = TRUE) %*% (Matrix::t(batch_mat) %*% grate[,"PTPRC"]))[,1]
      }
      if(pearson.family == "poisson"){
        yresp <- (yresp -  muhat) / sqrt(muhat )#+ muhat^2/10
      } else if (pearson.family == "nonzero_bernoulli"){
        nzprob <- 1-dpois(0, lambda = muhat)
        yresp <- ((yresp > 0) - nzprob) / (sqrt(nzprob * (1-nzprob)))
      }
      if(pearson.scale.max < Inf){
        yresp[yresp > pearson.scale.max*sd(yresp) + mean(yresp)] <- pearson.scale.max*sd(yresp) + mean(yresp)
      }
      yresp <- scale(yresp) 
    }
    if(verbose > 0){
      message(paste0(Sys.time(), ": begin fit"))
    }
    if(normalize_predictors_by_totalcounts){
      metagene_fit <-  
        spatial_durbin_nnls_wrap(y = yresp
                                 ,X = totalcount_norm(counts_matrix[,predictors], tc = totalcounts)
                                 ,Wrs = adjacency_matrix
                                 ,signs = markerslist[[jj]]$signs
                                 ,wts = prior_level_weights
                                 ,lag_response = lag_response
                                 ,lag_predictors = lag_predictors
                                 )
    } else {
      metagene_fit <-  
        spatial_durbin_nnls_wrap(y = yresp
                                 ,X = counts_matrix[,predictors]
                                 ,Wrs = adjacency_matrix
                                 ,signs = markerslist[[jj]]$signs
                                 ,wts = prior_level_weights
                                 ,lag_response = lag_response
                                 ,lag_predictors = lag_predictors
                                 )
      
    }
   
    if(verbose > 0){
      message(paste0(Sys.time(), ": finish fit"))
    } 
    metagene_fits[[ctclass]] <- list(
      bhat = metagene_fit$beta
      ,yhat = metagene_fit$xbeta[,1]
      ,y = yresp[,1]
    ) 
  }
  metagene_fits <- add_posteriormean(metagene_fits, v = "v1")
  return(metagene_fits)
}









