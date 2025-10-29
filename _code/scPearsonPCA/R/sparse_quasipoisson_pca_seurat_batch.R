#' Run sparse PCA Using Quasi-poisson pearson residuals computed by batch (skip dense normalization)
#'   
#' @description
#' 
#' Recreate Seurat RunPCA output from pearson residual normalization (described in Lause 2021), 
#' WITHOUT creating a dense cells x genes scaled data matrix.
#' This function uses Quasi-poisson instead of Negative Binomial distribution for computational convenience.
#' Computes PCA via SVD decomposition of genes x genes matrix instead of 
#' SVD decomposition of cells x genes dense matrix, which is useful for memory when cells x genes >> genes^2 .
#' 
#'
#' @param x A sparse genes x cells counts matrix 
#' @param obs a data.frame or data.table of metadata, containing the 'batch variable' and the 'cell id'.
#' @param batch_variable name of column in `obs` data frame  corresponding to 'batch'.
#' @param cellid_colname name of column in `obs` data frame corresponding to 'cell id'.
#' @param scale.max clip expression more than 'scale.max' standard deviations from the mean.  
#' Analogous to Seurat::ScaleData scale.max argument, with same default of "10".
#' @param reduction.key name of key for Seurat reduction object (if return_seurat_reduction=TRUE)
#' @param do.scale gene x gene cov matrix reflects as if pearson residuals are scaled by their standard deviation
#' @param do.center gene x gene cov matrix reflects as if the pearson residuals were centered by their mean
#' @param totalcounts  The total UMI counts across all gene targets in each cell.  This can be pre-computed by the user and passed in cases where the counts matrix `x` contains only a 
#' subset (of, say, highly variable) genes to be used in PCA. 
#' @param grate The gene frequencies.  This can be pre-computed by the user and passed in cases where the counts matrix `x` contains only a 
#' subset (of, say, highly variable) genes to be used in PCA. 
#' @param quasi_poisson_variance_inflation factor of overdispersion for Variance relative to the mean.  i.e., V(x) = E(x) * phi, where phi is the variance inflation factor
#' @param npcs the # of PCs
#' @param ncores the # of cores to use to project embeddings
#' @param verbose print out seurat-style PCA summary
#' @param block.size the # of cells to project into PCA space at a time. 
#' @param ndims.print the # of first PC's to print for verbose-seurat-style printout. 
#' @param nfeatures.print the # of top genes to show per PC for verbose-seurat-style printout. 
#' @param return_seurat_reduction if TRUE, returns a reduction.data object in Seurat format. 
#'
#'
#' @export  
sparse_quasipoisson_pca_seurat_batch <- function(x
                                                 ,obs
                                                 ,batch_variable
                                                 ,cellid_colname = "cell_ID"
                                                 ,scale.max=10
                                                 ,reduction.key="PC_"
                                                 ,do.scale=FALSE
                                                 ,do.center= FALSE
                                                 ,totalcounts = NULL
                                                 ,grate = NULL
                                                 ,quasi_poisson_variance_inflation = 1.01
                                                 ,npcs = 50
                                                 ,ncores=1
                                                 ,verbose=TRUE
                                                 ,block.size=10**4
                                                 ,ndims.print = 1:5
                                                 ,nfeatures.print = 30
                                                 ,return_seurat_reduction = TRUE
){
 
  x <- check_x_is_dgcmatrix(x)
  message(paste0(Sys.time(), ", computing PCA loadings..."))
  
  md <- data.table::copy(data.table::data.table(obs))
  if(!(cellid_colname) %in% names(md)){
    stopifnot(!is.null(rownames(obs)))
    cellid_colname <- "cell_ID"
    md[[cellid_colname]] <- rownames(obs)
  }
  if(cellid_colname!="cell_ID" & "cell_ID" %in% names(md)) md[["cell_ID"]] <- NULL
  setnames(md, old=cellid_colname, new="cell_ID")
  
  rm(obs); gc()
  md <- md[match(colnames(x),cell_ID)]
  
  md[,grpid__:=.GRP,by=c(batch_variable)]
  batch_mat <- Matrix::sparseMatrix(i = 1:nrow(md), j=md[["grpid__"]], x = 1
                                    ,dimnames = list(c(md[["cell_ID"]])
                                                     ,c(md[,head(.SD, 1),by=grpid__][order(grpid__)][[batch_variable]])
                                    )
  )
  
  phi <- check_phi(x, quasi_poisson_variance_inflation)
  nn <- ncol(x)
  totalcounts <- check_totalcounts(x, totalcounts)
  
  ### expression rates by batch (\hat{p})
  grate <- check_grate_batch(x, grate, batch_mat)
  
  if (.Platform$OS.type == "windows" && ncores > 1L) {
    warning("mclapply() runs serially on Windows; mc.cores ignored.")
    ncores <- 1L
  } 
  ###  
#  cellnames_vec <- findInterval(seq(x@x)-1, x@p[-1]) + 1
  cellnames_vec <- rep.int(seq_len(ncol(x)), diff(x@p))
  batch_subset <- vector(mode='list',length=ncol(grate)) 
  for(ii in 1:ncol(grate)){
    batch_idx <- md[grpid__==ii,which=TRUE]                    ## which cell_IDs belong to batch 'ii'
    batch_subset[[ii]] <- which(cellnames_vec %in% batch_idx)  ## which elements of x matrix belong to batch 'ii'
  } 
  
  
  ytilde <- x ## a copy
  for(ii in 1:ncol(grate)){
    message(ii)
    grate_inp <- sqrt(1/(grate[,ii]*phi))
    grate_inp[is.infinite(grate_inp)] <- 0 ## avoid Inf * 0 / NA error
    names(grate_inp) <- rownames(grate)
    grate_subset <- rownames(x)[(x@i + 1)][batch_subset[[ii]]] ## genes corresponding to batch 'ii'
    ytilde@x[batch_subset[[ii]]] <- ytilde@x[batch_subset[[ii]]] * grate_inp[grate_subset]
    gc()
  }
  
  ### diagonal matrix (1/totalcounts)
  root_tc_diag <- Matrix::Diagonal(x = sqrt(1/totalcounts))
  dimnames(root_tc_diag) <- list(colnames(x) , colnames(x))
  
  ytilde <- ytilde %*% root_tc_diag
  
  ###  
  root_tc_col <- 
    Matrix::Diagonal(x=sqrt(totalcounts)) %*% batch_mat
  dimnames(root_tc_col) <- list(colnames(x) , colnames(batch_mat))
  
  root_grate_phi_row <- grate
  root_grate_phi_row@x <- sqrt(root_grate_phi_row@x / phi)
  
  ytilde_muhat <- ytilde %*% root_tc_col
  ytilde_muhat <- ytilde_muhat %*%  Matrix::t(root_grate_phi_row)
  
  mean_sqrt_tc_batch <- 
    Matrix::Diagonal(x=Matrix::colSums(Matrix::Diagonal(x=sqrt(totalcounts)) %*% batch_mat)) %*%
    Matrix::Diagonal(x=1/Matrix::colSums(batch_mat))
  dimnames(mean_sqrt_tc_batch) <- list(colnames(batch_mat) , colnames(batch_mat))
   
   
  batch_rate <- Matrix::colSums(batch_mat)
  batch_rate <- batch_rate / sum(batch_rate)
  
  mean.pearson_residual <- (root_grate_phi_row %*% mean_sqrt_tc_batch %*% batch_rate)[,1]
  mean.pearson_residual <- Matrix::rowMeans(ytilde) - mean.pearson_residual
  
  inner_tc <-  Matrix::Diagonal(x = sqrt(totalcounts)
                                ,names =colnames(x)) %*% batch_mat
  dimnames(inner_tc) <- list(colnames(x), colnames(batch_mat))
  inner_tc <- Matrix::crossprod(inner_tc)
  
  #### The cross product of pearson residuals:
  #### [y/sqrt(v(y)) -  muhat/sqrt(v(y))] [y/sqrt(v(y)) -  muhat/sqrt(v(y))]^T
  qp <- ytilde %*% Matrix::t(ytilde) - 
    ytilde_muhat - Matrix::t(ytilde_muhat) + 
    root_grate_phi_row %*% inner_tc %*% Matrix::t(root_grate_phi_row)
  
  gc()
  
  #### Note: diagonal elements of the crossproduct matrix 'qp' give \sum y_i^2 for each gene             
  #### This can be used to get the standard deviation of pearson residuals
  #### i.e., V(x) = E(x^2) - E(x)^2
  sd.pearson_residual <- sqrt((Matrix::diag(qp) - 
                                 nn * mean.pearson_residual^2 ) / 
                                (nn-1))
  
  if(scale.max < Inf){
    message_parallel(paste0(Sys.time(), ", computing clipped values of pearson residuals > scale.max=", scale.max, " sd's above the mean pearson residual")) 
    clip_vals <- scale.max * sd.pearson_residual + mean.pearson_residual
    xvec <- ytilde@x   
    names(xvec) <- rownames(ytilde)[(ytilde@i + 1)] 
    max_val_vec <- clip_vals[names(xvec)] 
    
    for(ii in 1:ncol(grate)){
      message(ii)
      grate_subset <- rownames(x)[(x@i + 1)][batch_subset[[ii]]] ## genes corresponding to batch 'ii'
      max_val_vec[batch_subset[[ii]]] <- max_val_vec[batch_subset[[ii]]] + 
        sqrt(grate[grate_subset,ii]/phi) * sqrt(totalcounts[cellnames_vec[batch_subset[[ii]]]])
      gc()
    }
    xvec <- pmin(xvec, max_val_vec) 
    ytilde@x <- xvec
    
    ytilde_muhat <- ytilde %*% root_tc_col
    ytilde_muhat <- ytilde_muhat %*%  Matrix::t(root_grate_phi_row)
    
    #### The cross product of pearson residuals:
    #### [y/sqrt(v(y)) -  muhat/sqrt(v(y))] [y/sqrt(v(y)) -  muhat/sqrt(v(y))]^T
    qp <- ytilde %*% Matrix::t(ytilde) - 
      ytilde_muhat - Matrix::t(ytilde_muhat) + 
      root_grate_phi_row %*% inner_tc %*% Matrix::t(root_grate_phi_row)
    
  }
  
  if(do.center){
    
    message_parallel(paste0(Sys.time(), ", centering pearson residuals")) 
    mean.pearson_residual.clipped <- (root_grate_phi_row %*% mean_sqrt_tc_batch %*% batch_rate)[,1]
    mean.pearson_residual.clipped <- Matrix::rowMeans(ytilde) - mean.pearson_residual.clipped
    
    inner_prod <- nn*mean.pearson_residual.clipped %*% t(mean.pearson_residual)  
    
    qp <-  
      (qp - inner_prod - Matrix::t(inner_prod) + 
         nn*Matrix::tcrossprod(mean.pearson_residual))
    
    gc()
  }
  
  if(do.scale){
    message_parallel(paste0(Sys.time(), ", scaling pearson residuals")) 
    qp <- Matrix::Diagonal(x = 1/sd.pearson_residual) %*% qp %*% Matrix::Diagonal(x = 1/sd.pearson_residual)
    dimnames(qp) <- list(names(sd.pearson_residual), names(sd.pearson_residual))
  } 
  
  #### Get the eigenvectors / loadings 
  svdd <- RSpectra::svds(qp, k = min(nrow(qp), npcs))
  feature.loadings <- svdd$u
  rownames(feature.loadings) <- rownames(x)
  colnames(feature.loadings) <- paste0(reduction.key, 1:npcs)
  
  message_parallel(paste0(Sys.time(), ", computing PCA cell embeddings using ", ncores, " threads.")) 
  
  gc() 
  ### Need to get the cell embeddings.
  ### scale one piece at a time and project 
  splitdt <- data.table(i=1:ncol(x))[,ss:=floor(i/block.size) + 1L]
  splitidx <- lapply(split(splitdt, by="ss"), function(xx){ xx[["i"]]})
  cell.embeddings <-  
    parallel::mclapply(1:length(splitidx)
                       ,function(kk
                                 ,doscale=do.scale
                                 ,docenter=do.center
                                 ,centers = mean.pearson_residual
                                 ,scales = sd.pearson_residual
                       ){
                         scx <-  
                           Matrix::t(ytilde[,splitidx[[kk]],drop=FALSE]) - 
                           Matrix::Diagonal(x=sqrt(totalcounts[splitidx[[kk]]])) %*% batch_mat[splitidx[[kk]],]  %*%
                           Matrix::t(root_grate_phi_row)
                         
                         if(doscale & docenter){
                           scx <- scale(as.matrix(scx), center = centers, scale = scales) 
                         } else if (docenter) {
                           scx <- scale(as.matrix(scx), center = centers, scale = FALSE) 
                           
                         } else if (doscale) {
                           scx <- scale(as.matrix(scx), center = FALSE, scale = scales) 
                         }
                         return(scx %*% feature.loadings)
                       },mc.cores=ncores
    ) 
  cell.embeddings <- do.call(rbind, cell.embeddings)
  
  if(return_seurat_reduction){
    sdev <- sqrt(svdd$d / (ncol(x)-1)) 
    reduction.data <- 
      Seurat::CreateDimReducObject(
        embeddings = cell.embeddings
        ,loadings = feature.loadings
        ,assay = "RNA"
        ,stdev = sdev
        ,key = reduction.key
      )
    
    if (verbose) {
      msg <- utils::capture.output(print(
        x = reduction.data,
        dims = ndims.print,
        nfeatures = nfeatures.print
      ))
      message(paste(msg, collapse = '\n'))
    }
    
    return(list(reduction.data = reduction.data
                ,mean.pearson_residual = mean.pearson_residual
                ,sd.pearson_residual = sd.pearson_residual
                ,scale.max = scale.max
                ,do.scale = do.scale
                ,do.center = do.center
                ,grate = grate
                ,phi = phi))
      
  } else {
    return(list(cell.embeddings = cell.embeddings
                ,feature.loadings = feature.loadings
                ,mean.pearson_residual = mean.pearson_residual
                ,sd.pearson_residual = sd.pearson_residual
                ,scale.max = scale.max
                ,do.scale = do.scale
                ,do.center = do.center
                ,grate = grate
                ,phi = phi))
    
  }
}



