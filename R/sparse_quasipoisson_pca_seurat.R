#' Run sparse PCA Using Quasi-poisson pearson residuals (skip dense normalization)
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
#' @param totalcounts  The total UMI counts across all gene targets in each cell.  This can be pre-computed by the user and passed in cases where the counts matrix `x` contains only a 
#' subset (of, say, highly variable) genes to be used in PCA. 
#' @param grate The gene frequencies.  This can be pre-computed by the user and passed in cases where the counts matrix `x` contains only a 
#' subset (of, say, highly variable) genes to be used in PCA. 
#' @param scale.max clip expression more than 'scale.max' standard deviations from the mean.  
#' Analogous to Seurat::ScaleData scale.max argument, with same default of "10".
#' @param reduction.key Seurat key
#' @param do.scale gene x gene cov matrix reflects as if pearson residuals are scaled by their standard deviation
#' @param do.center gene x gene cov matrix reflects as if the pearson residuals were centered by their mean
#' @param quasi_poisson_variance_inflation factor of overdispersion for Variance relative to the mean.  i.e., V(x) = E(x) * phi, where phi is the variance inflation factor
#' @param npcs the # of PCs
#' @param ncores the # of cores to use to project embeddings
#' @param verbose print out seurat-style PCA summary
#' @param block.size the # of cells to project into PCA space at a time. Can reduce block.size for lower-memory computing.
#' @param ndims.print the # of first PC's to print for verbose-seurat-style printout. 
#' @param nfeatures.print the # of top genes to show per PC for verbose-seurat-style printout. 
#' @param return_seurat_reduction if TRUE, returns a reduction.data object in Seurat format. 
#'
#'
#' @export  
sparse_quasipoisson_pca_seurat <- function(x
                                           ,totalcounts = NULL
                                           ,grate = NULL
                                           ,scale.max=10
                                           ,reduction.key="PC_"
                                           ,do.scale=TRUE
                                           ,do.center= TRUE
                                           ,quasi_poisson_variance_inflation = 1.01
                                           ,npcs = 50
                                           ,ncores=1
                                           ,verbose=TRUE
                                           ,block.size=10**4
                                           ,ndims.print = 1:5
                                           ,nfeatures.print = 30
                                           ,return_seurat_reduction = TRUE
                                           ,estimate_dispersion_from_poisson = FALSE
                                           ,only_return_sds = FALSE
){
  
  x <- check_x_is_dgcmatrix(x)
  nn <- ncol(x)
  phi <- quasi_poisson_variance_inflation   
  totalcounts <- check_totalcounts(x, totalcounts)
  grate <- check_grate(x, grate)
  
  message_parallel(paste0(Sys.time(), ", computing PCA loadings..."))
  ### diagonal matrix (1/(phi*estimated gene frequency))
  root_grate_phi_diag <- Matrix::Diagonal(x = sqrt(1/(grate * phi))
                                          ,names =names(grate))
  
  ### diagonal matrix (1/totalcounts)
  root_tc_diag <- Matrix::Diagonal(x = sqrt(1/totalcounts)
                                   ,names =colnames(x))
  
  ### y / sqrt(v(y)) 
  ytilde <- root_grate_phi_diag %*% x %*% root_tc_diag ## sparse  n x g
  
  ### gene-wise average of y/sqrt(v(y)) - muhat/sqrt(v(y))
  ### Note: muhat / sqrt(v(y)) = muhat/sqrt(muhat*phi) = sqrt(muhat/phi)
  mean.pearson_residual <- Matrix::rowMeans(ytilde) - 
    mean(sqrt(totalcounts)) * sqrt(grate/phi)
  
  
  root_grate_phi_row <- Matrix::Matrix(data = sqrt(grate / phi)
                                       ,nrow=1
                                       ,dimnames =list(c(), names(grate)))
  root_tc_col <- 
    Matrix::Matrix(sqrt(totalcounts), ncol = 1
                   ,dimnames = list(c(colnames(x)), c()))
  
  ytilde_muhat <- ytilde %*% root_tc_col ### g x 1
  ytilde_muhat <- ytilde_muhat %*% root_grate_phi_row ### g x g
  
  
  message_parallel(paste0(Sys.time(), ", computing cross-product of pearson residuals.")) 
  #### The cross product of pearson residuals:
  #### [y/sqrt(v(y)) -  muhat/sqrt(v(y))] [y/sqrt(v(y)) -  muhat/sqrt(v(y))]^T
  qp <- ytilde %*% Matrix::t(ytilde) - 
    ytilde_muhat - Matrix::t(ytilde_muhat) + 
    sum(root_tc_col^2) * Matrix::t(root_grate_phi_row) %*% root_grate_phi_row
  
  if(estimate_dispersion_from_poisson & phi==1){
    ##browser()
    phi <- Matrix::diag(qp) / (nn - 1) ## get glm-like estimate of overdispersion factor phi
    ### diagonal matrix (1/(phi*estimated gene frequency))
    root_grate_phi_diag <- Matrix::Diagonal(x = sqrt(1/(grate * phi))
                                            ,names =names(grate))
    
    ### diagonal matrix (1/totalcounts)
    root_tc_diag <- Matrix::Diagonal(x = sqrt(1/totalcounts)
                                     ,names =colnames(x))
    
    ### y / sqrt(v(y)) 
    ytilde <- root_grate_phi_diag %*% x %*% root_tc_diag ## sparse  n x g
    
    ### gene-wise average of y/sqrt(v(y)) - muhat/sqrt(v(y))
    ### Note: muhat / sqrt(v(y)) = muhat/sqrt(muhat*phi) = sqrt(muhat/phi)
    mean.pearson_residual <- Matrix::rowMeans(ytilde) - 
      mean(sqrt(totalcounts)) * sqrt(grate/phi)
    
    
    root_grate_phi_row <- Matrix::Matrix(data = sqrt(grate / phi)
                                         ,nrow=1
                                         ,dimnames =list(c(), names(grate)))
    root_tc_col <- 
      Matrix::Matrix(sqrt(totalcounts), ncol = 1
                     ,dimnames = list(c(colnames(x)), c()))
    
    ytilde_muhat <- ytilde %*% root_tc_col ### g x 1
    ytilde_muhat <- ytilde_muhat %*% root_grate_phi_row ### g x g
    
    
    message_parallel(paste0(Sys.time(), ", computing cross-product of pearson residuals.")) 
    #### The cross product of pearson residuals:
    #### [y/sqrt(v(y)) -  muhat/sqrt(v(y))] [y/sqrt(v(y)) -  muhat/sqrt(v(y))]^T
    qp <- ytilde %*% Matrix::t(ytilde) - 
      ytilde_muhat - Matrix::t(ytilde_muhat) + 
      sum(root_tc_col^2) * Matrix::t(root_grate_phi_row) %*% root_grate_phi_row
      
  } 
  #### Note: diagonal elements of the crossproduct matrix 'qp' give \sum y_i^2 for each gene             
  #### This can be used to get the standard deviation of pearson residuals
  #### i.e., V(x) = E(x^2) - E(x)^2
  sd.pearson_residual <- sqrt((Matrix::diag(qp) - 
                                 nn * mean.pearson_residual^2 ) / 
                                (nn-1))
 
  if(only_return_sds){
    return(list(sd.pearson_residual = sd.pearson_residual))
  } 
  #### We can infer the cell/gene specific values 
  #### which are scale.max SD's above the mean and clip them 
  if(scale.max < Inf){
    message_parallel(paste0(Sys.time(), ", computing clipped values of pearson residuals > scale.max=", scale.max, " sd's above the mean pearson residual")) 
    clip_vals <- scale.max * sd.pearson_residual + mean.pearson_residual
    xvec <- ytilde@x   
    names(xvec) <- rownames(ytilde)[(ytilde@i + 1)] 
    max_val_vec <- clip_vals[names(xvec)] 
    cellnames_max_val_vec <- colnames(ytilde)[findInterval(seq(ytilde@x)-1, ytilde@p[-1]) + 1]
    
    ## pearson residual = y' - mu'
    ## y' = y/sqrt(p_g n_c phi ); mu' = sqrt(p_g n_c/phi)
    ## clipping: y' - mu' > c ==> y' > mu' + c
    ## clip values of y' > mu' + c, to mu' + c
    if(length(phi)==1){
      max_val_vec <- max_val_vec + 
        sqrt(grate[names(max_val_vec)]/phi) * sqrt(totalcounts[cellnames_max_val_vec])
        
    } else if(length(phi)==length(grate)){
      phi_val_vec <- phi[names(xvec)]  
      max_val_vec <- max_val_vec + 
        sqrt(grate[names(max_val_vec)]/phi_val_vec[names(xvec)]) * sqrt(totalcounts[cellnames_max_val_vec])
    } else {
      stop(paste0("overdispersion factor has different length (", length(phi), ") than the number of genes"))
    }
    xvec <- pmin(xvec, max_val_vec) 
    
    ytilde@x <- xvec
    rm(list=c("xvec", "max_val_vec", "cellnames_max_val_vec")); gc()
    
    ytilde_muhat <- ytilde %*% root_tc_col ### g x 1
    ytilde_muhat <- ytilde_muhat %*% root_grate_phi_row ### g x g
    
    qp <- ytilde %*% Matrix::t(ytilde) - 
      ytilde_muhat - Matrix::t(ytilde_muhat) + 
      sum(root_tc_col^2) * Matrix::t(root_grate_phi_row) %*% root_grate_phi_row
  }
  
  if(do.center){
    
    message_parallel(paste0(Sys.time(), ", centering pearson residuals")) 
    mean.pearson_residual.clipped  <- 
      Matrix::rowMeans(ytilde) - 
      mean(sqrt(totalcounts)) * sqrt(grate/phi)
    
    inner_prod <- nn*mean.pearson_residual.clipped %*% t(mean.pearson_residual)  
    
    qp <-  
      (qp - inner_prod - Matrix::t(inner_prod) + 
         nn*Matrix::tcrossprod(mean.pearson_residual))
    
  }
  
  if(do.scale){
    message_parallel(paste0(Sys.time(), ", scaling pearson residuals")) 
 #   qp <- diag(1/sd.pearson_residual) %*% qp %*% diag(1/sd.pearson_residual)
    qp <- Matrix::Diagonal(x = 1/sd.pearson_residual, names = names(sd.pearson_residual)) %*% qp %*% Matrix::Diagonal(x = 1/sd.pearson_residual, names = names(sd.pearson_residual))
  } 
  
  #### Get the eigenvectors / loadings 
  svdd <- RSpectra::svds(qp, k = npcs)
  feature.loadings <- svdd$u
  rownames(feature.loadings) <- rownames(x)
  colnames(feature.loadings) <- paste0(reduction.key, 1:npcs)
  
  message_parallel(paste0(Sys.time(), ", computing PCA cell embeddings using ", ncores, " threads.")) 
  
  gc() 
  ### Need to get the cell embeddings.
  ### scale one piece at a time and project 
  splitdt <- data.table::data.table(i=1:ncol(x))[,ss:=floor(i/block.size) + 1L]
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
                           matrix(sqrt(totalcounts[splitidx[[kk]]]),ncol=1) %*% matrix(sqrt(grate/phi),nrow=1)
                         
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
      
  }  else {
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
