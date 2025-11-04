




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
#' @param training_ids optional vector of cell IDs to use for training NNLS metagene models.
#' @param modeltype Right now, just 'nnls' (non-negative least squares) is supported.
#' @param lag_response Should the (weighted) average of the index gene across similar cells be included as a predictor in the regression model? default is TRUE
#' @param lag_predictors Should the (weighted) averages of predictor genes across similar cells be included as predictors in the regression model? default is TRUE
#' @param pearson.response Should the response variable be the pearson residuals of the index variable?
#' @param pearson.scale.max If `pearson.response = TRUE`, have large pearson residuals > `pearson.scale.max` SDs above the mean truncated to this amount.
#' @param pearson.family "poisson" pearson residuals (recommended), or "nonzero_bernoulli" (experimental).
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
                         ,training_ids = NULL
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
                         ,...
  ){
 
  stopifnot(inherits(markerslist, "markerslist")) 
  
  # 
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

  rnapreds <- grep("_protein$", colnames(counts_matrix), value=TRUE, invert=TRUE)
  proteinpreds <- setdiff(colnames(counts_matrix), rnapreds)
  
  ### get totalcounts if not specified 
  if(length(rnapreds) > 0){
    if(is.null(totalcounts)){
      totalcounts <- Matrix::rowSums(counts_matrix[,rnapreds,drop=FALSE]) 
    }
    stopifnot(length(totalcounts) == nrow(counts_matrix))
    if(!is.null(names(totalcounts)) && !is.null(rownames(counts_matrix))){
      stopifnot("ids of 'totalcounts' dont match ids in 'counts_matrix'" = 
                  all(names(totalcounts) %in% rownames(counts_matrix))) 
      totalcounts <- totalcounts[rownames(counts_matrix)] 
    }
  }
  
  ### check length and names of prior_level_weights 
  if(!is.null(prior_level_weights)){
    stopifnot(length(prior_level_weights) == nrow(counts_matrix))
    if(!is.null(names(prior_level_weights)) && !is.null(rownames(counts_matrix))){
      stopifnot("ids of 'prior_level_weights' dont match ids in 'counts_matrix'" = 
                  all(names(prior_level_weights) %in% rownames(counts_matrix))) 
      prior_level_weights <- prior_level_weights[rownames(counts_matrix)] 
    } else {
      names(prior_level_weights) <- rownames(counts_matrix)
    }
  } 
 
  ### compute frequency of each gene target (for use in pearson residuals calculation) 
  if(!is.null(batch_variable)){
    ## make batch_mat
    md <- data.table::copy(data.table::data.table(obs))
    if(!(cellid_colname) %in% names(md)){
      stopifnot(!is.null(rownames(obs)))
      cellid_colname <- "cell_ID"
      md[[cellid_colname]] <- rownames(obs)
    }
    if(cellid_colname!="cell_ID" & "cell_ID" %in% names(md)) md[["cell_ID"]] <- NULL
    data.table::setnames(md, old=cellid_colname, new="cell_ID")
    
    rm(obs); gc()
    md <- md[match(rownames(counts_matrix),cell_ID)]
    
    md[,grpid__:=.GRP,by=c(batch_variable)]
    batch_mat <- Matrix::sparseMatrix(j = 1:nrow(md), i=md[["grpid__"]], x = 1
                                      ,dimnames = list(c(md[,head(.SD, 1),by=grpid__][order(grpid__)][[batch_variable]])
                                                       ,c(md[["cell_ID"]]))
                                                       
                                      )
    ## make batch_mat
    if(is.null(gene_wise_frequency)){
      ### expression rates by batch (\hat{p})
      grate <- Matrix::Diagonal(x=1/Matrix::rowSums(batch_mat), names = rownames(batch_mat)) %*% batch_mat %*%  counts_matrix 
      grate <- Matrix::Diagonal(x = 1/Matrix::rowSums(grate),names=rownames(batch_mat)) %*% grate 
    } else {
      if(all(rownames(batch_mat) %in% colnames(gene_wise_frequency))){
        grate <- Matrix::t(gene_wise_frequency[,rownames(batch_mat)])
      } else {
        if(!all(colnames(batch_mat)) %in% rownames(gene_wise_frequency)){
          grate <- gene_wise_frequency[colnames(batch_mat),]
        }
      }
    } 
  }  else {
    if(is.null(gene_wise_frequency)){
      gene_wise_frequency <- Matrix::colSums(counts_matrix)
      gene_wise_frequency <- gene_wise_frequency / sum(gene_wise_frequency)
    }
  }
  
 
  if(length(allgenes) < ncol(counts_matrix)){
    counts_matrix <- counts_matrix[,allgenes]  
    gc()
  } 
  rnapreds <- intersect(rnapreds, allgenes)
  proteinpreds <- intersect(proteinpreds, allgenes)
  
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

  #### Pre-compute some matrix multiplications that will be 
  #### re-used across celltype classes
  Xall <- counts_matrix
  if(normalize_predictors_by_totalcounts){
    message(paste0("totalcount norm (all): ", Sys.time()))
    if(length(rnapreds) > 0){
      Xallrna <- totalcount_norm(Xall[,c(rnapreds),drop=FALSE], tc = totalcounts)
      if(length(proteinpreds) > 0){
        Xall <- cbind(Xallrna, Xall[,proteinpreds,drop=FALSE]) 
      } else {
        Xall <- Xallrna
      }
      rm(Xallrna); gc()
    }
  }
  if(!is.null(adjacency_matrix)){
    if(lag_predictors){
      message(paste0("lagging predictors (all): ", Sys.time()))
      WXall <- adjacency_matrix %*% Xall
    }
  }
 
  #### Loop over cell type classes to get NNLS fits and scores 
  for(jj in 1:length(markerslist)){
    ctclass <- names(markerslist)[jj]
    message(paste0("Modeling ", ctclass, ". ",  Sys.time()))
    index_markers <- markerslist[[ctclass]][["index_marker"]] 
    predictors <- markerslist[[ctclass]]$predictors
    signs <- rep(1, length(predictors)) 
    if(markerslist[[ctclass]]$use_offclass_markers_as_negative_predictors){
       offclass <- unname(unlist(lapply(markerslist[-jj], "[[", "predictors")))
       offclass <- setdiff(offclass, predictors)
       offclass_signs <- rep(-1, length(offclass))
       predictors <- c(predictors, offclass)
       signs <- c(signs, offclass_signs)
    }
    yrespmat <- vector(mode = 'list',length=length(index_markers)) 
    
    for(index_marker in index_markers){
      yresp <- counts_matrix[,index_marker]
      if(!grepl("_protein$", index_marker)){
        if(pearson.response){
          ### summary statistics for normalization and 
          ### (optional) pearson standardization of response variable used for fitting metagene
          if(is.null(batch_variable)){
              muhat <- gene_wise_frequency[index_marker] * totalcounts
          } else {
            muhat <- (Matrix::Diagonal(x = totalcounts, names = TRUE) %*% (Matrix::t(batch_mat) %*% grate[,index_marker]))[,1]
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
        }
      } 
      yrespmat[[index_marker]] <- scale(yresp)
    }
    yrespmat <- do.call(cbind, yrespmat)
    colnames(yrespmat) <- index_markers
   
    Xsigned <- Xall[,predictors] %*%  Matrix::Diagonal(x = signs, names = predictors)
    
    if(!is.null(adjacency_matrix)){
      if(lag_predictors){
        message(paste0("lagging predictors: ", Sys.time()))
        #WXmat <- adjacency_matrix %*% Xsigned
        #WXmat <- adjacency_matrix %*% Xsigned
        WXmat <- WXall[,predictors] %*% Matrix::Diagonal(x = signs, names = predictors)
        colnames(WXmat) <- paste0("lag.", colnames(WXmat))
        Xsigned <- cbind(Xsigned, WXmat)
        signs <- c(signs, signs)
        rm(WXmat); gc()
      }
    }
   
    if(lag_response){
      message(paste0("lagging response: ", Sys.time()))
      
      Wy <- adjacency_matrix %*% yrespmat
      colnames(Wy) <- paste0("lag.y.", colnames(yrespmat))
      Xsigned <- cbind(Wy, Xsigned)
      signs <- c(rep(1,ncol(Wy)), signs)
      rm(Wy)
    }  
    
    message(paste0("cross-prod calcualtions: ", Sys.time()))
    
    if(!is.null(training_ids)){
      if(!is.null(prior_level_weights)){
        a <- Matrix::crossprod(Matrix::Diagonal(x=sqrt(prior_level_weights[training_ids])) %*% Xsigned[training_ids,])
        b <- Matrix::crossprod(Matrix::Diagonal(x=prior_level_weights[training_ids]) %*% Xsigned[training_ids] , yrespmat[training_ids,,drop=FALSE])
      } else {
        a <- Matrix::crossprod(Xsigned[training_ids,])
        b <- Matrix::crossprod(Xsigned[training_ids,], yrespmat[training_ids,,drop=FALSE])
      }
    } else {
      if(!is.null(prior_level_weights)){
        a <- Matrix::crossprod(Matrix::Diagonal(x=sqrt(prior_level_weights)) %*% Xsigned)
        b <- Matrix::crossprod(Matrix::Diagonal(x=prior_level_weights) %*% Xsigned , yrespmat)
      } else {
        a <- Matrix::crossprod(Xsigned)
        b <- Matrix::crossprod(Xsigned, yrespmat)
      }
    }
    gc()
    unconstrained <- FALSE 
    for(index_marker in colnames(yrespmat)){
      message(paste0("fitting model to index marker: ", index_marker, ". ", Sys.time()))
      rm_these_predictors <- which(colnames(a) %in% paste0(c("", "lag."),index_marker))
      rm_these_predictors <- sort(c(rm_these_predictors, which(colnames(a) %in% paste0(c("lag."),setdiff(index_markers, index_marker)))))
      #colnames(a)[rm_these_predictors]
      if(unconstrained){
        beta <- solve(as.matrix(a[-c(rm_these_predictors),-c(rm_these_predictors),drop=FALSE])
                      ,as.matrix(b[-c(rm_these_predictors),index_marker,drop=FALSE]))
      } else {
        beta <- RcppML::nnls(as.matrix(a[-c(rm_these_predictors),-c(rm_these_predictors),drop=FALSE])
                             ,as.matrix(b[-c(rm_these_predictors),index_marker,drop=FALSE]))
      }
      xbeta <- Xsigned[,-c(rm_these_predictors),drop=FALSE] %*% beta ## get prediction
      rownames(beta) <- colnames(Xsigned[,-c(rm_these_predictors),drop=FALSE])
     
      beta <- beta * signs[-c(rm_these_predictors)] ## convert back
      
      metagene_fits[[ctclass]][[index_marker]] <- list(
        bhat = beta
        ,yhat = xbeta[,1]
        ,y = yrespmat[,index_marker]
      )
    } 
  }
  gc() 
  ### Summarize the posterior mean scores for each index marker in each celltype class
  metagene_fits <- add_posteriormean(metagene_fits, wts = prior_level_weights, v = "v1")
  
  ### For multiple index_markers, use the non-negative sparse-pca embedding to 
  ### summarize a metagene score for each cell type class in one dimension
  ### Note - Non-negative sparse PCA (used to summarize scores in 1 dimension across index gene)
  ###        uses an em-algorithm for loadings.  Need to set a seed for full reproducibility
  orig_seed <- .Random.seed
  on.exit({.Random.seed <<- orig_seed})
  nsprcomp_seed <- 123
  set.seed(nsprcomp_seed) 
  nnpc_for_multi_index <- TRUE
  if(nnpc_for_multi_index){
    for(ct in 1:length(metagene_fits)){
      set.seed(nsprcomp_seed)
      yct <- do.call(cbind, lapply(metagene_fits[[ct]], "[[", "ypost"))
      nspc <- nsprcomp::nsprcomp(yct, scale.=TRUE, nneg = TRUE, ncomp = 1)
      ypostct <- (scale(yct) %*% nspc$rotation[,1])[,1]
      
      yct <- do.call(cbind, lapply(metagene_fits[[ct]], "[[", "yhat"))
      nspc <- nsprcomp::nsprcomp(yct, scale.=TRUE, nneg = TRUE, ncomp = 1)
      yhatct <- (scale(yct) %*% nspc$rotation[,1])[,1]
      
      yct <- do.call(cbind, lapply(metagene_fits[[ct]], "[[", "y"))
      nspc <- nsprcomp::nsprcomp(yct, scale.=TRUE, nneg = TRUE, ncomp = 1)
      yobs <- (scale(yct) %*% nspc$rotation[,1])[,1]
    
    
      bhat <- do.call(cbind, lapply(metagene_fits[[ct]], "[[", "bhat")) 
      bhatcomb <- lapply(metagene_fits[[ct]], "[[", "bhat")
      bhatcomb <- lapply(names(bhatcomb), function(ii){
        ret <- data.table::data.table(bhatcomb[[ii]], keep.rownames = TRUE)
        ret[rn=="lag.y",rn:=paste0("lag.",ii)][,index_marker:=ii]
        return(ret)
      })
      bhatcomb <- data.table::dcast(data.table::rbindlist(bhatcomb), rn ~ index_marker, value.var="V1",fill=0)
      rn <- bhatcomb[[1]]
      bhatcomb <- as(bhatcomb[,-1], "sparseMatrix")
      rownames(bhatcomb) <- rn
      bhatct <- bhatcomb %*% Matrix::Diagonal(x = apply(yct, 2, sd) , names = TRUE) %*% nspc$rotation[,1]
      bhatct <- bhat[order(-bhat[,1]),,drop=FALSE]
      metagene_fits[[ct]][["yhat"]] <- yhatct
      metagene_fits[[ct]][["y"]] <- yobs
      metagene_fits[[ct]][["ypost"]] <- ypostct
      metagene_fits[[ct]][["bhat"]] <- bhatct
      metagene_fits[[ct]][["nnpcfit"]] <- nspc
    }
  }
  return(metagene_fits)
}









