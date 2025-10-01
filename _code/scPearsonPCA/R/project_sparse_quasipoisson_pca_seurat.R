
#' Project PCA onto new data 
#' 
#' @description
#' 
#' Take the results generated from `sparse_quasipoisson_pca_seurat` and project them onto new cells.
#' 
#'
#' @param x A sparse genes x cells counts matrix 
#' @param pcaobj generated from `sparse_quasipoisson_pca_seurat`  
#' Analogous to Seurat::ScaleData scale.max argument, with same default of "10".
#' @param ncores the # of cores to use to project embeddings
#' @param verbose print out seurat-style PCA summary
#' @param block.size the # of cells to project into PCA space at a time. Can reduce block.size for lower-memory computing.
#' @param ndims.print the # of first PC's to print for verbose-seurat-style printout. 
#' @param nfeatures.print the # of top genes to show per PC for verbose-seurat-style printout. 
#' @param return_seurat_reduction if TRUE, returns a reduction.data object in Seurat format. 
#'
#'
#' @export  
project_sparse_quasipoisson_pca_seurat <- function(x
                                                   ,pcaobj
                                                   ,ncores=1
                                                   ,verbose=TRUE
                                                   ,block.size=10**4
                                                   ,ndims.print = 1:5
                                                   ,nfeatures.print = 30
                                                   ,return_seurat_reduction = TRUE
){
  
  nn <- ncol(x)
  phi <- pcaobj$phi# quasi_poisson_variance_inflation   
  totalcounts <- Matrix::colSums(x)
  grate <- pcaobj$grate #Matrix::rowSums(x) 
  
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
  mean.pearson_residual <- pcaobj$mean.pearson_residual
  
  
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
  
  
  #### Note: diagonal elements of the crossproduct matrix 'qp' give \sum y_i^2 for each gene             
  #### This can be used to get the standard deviation of pearson residuals
  #### i.e., V(x) = E(x^2) - E(x)^2
  sd.pearson_residual <- pcaobj$sd.pearson_residual
  
  #### We can infer the cell/gene specific values 
  #### which are scale.max SD's above the mean and clip them 
  scale.max <- pcaobj$scale.max
  do.center <- pcaobj$do.center
  do.scale <- pcaobj$do.scale
  if(scale.max < Inf){
    message_parallel(paste0(Sys.time(), ", computing clipped values of pearson residuals > scale.max=", scale.max, " sd's above the mean pearson residual")) 
    clip_vals <- scale.max * sd.pearson_residual + mean.pearson_residual
    xvec <- ytilde@x   
    names(xvec) <- rownames(ytilde)[(ytilde@i + 1)] 
    max_val_vec <- clip_vals[names(xvec)] 
    cellnames_max_val_vec <- colnames(ytilde)[findInterval(seq(ytilde@x)-1, ytilde@p[-1]) + 1]
    
    max_val_vec <- max_val_vec + 
      sqrt(grate[names(max_val_vec)]/phi) * sqrt(totalcounts[cellnames_max_val_vec])
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
    qp <- diag(1/sd.pearson_residual) %*% qp %*% diag(1/sd.pearson_residual)
  } 
  
  #### Get the eigenvectors / loadings 
  feature.loadings <- pcaobj$reduction.data@feature.loadings
  sdev <- pcaobj$reduction.data@stdev
  reduction.key <- pcaobj$reduction.data@key
  
  message_parallel(paste0(Sys.time(), ", computing PCA cell embeddings using ", ncores, " threads.")) 
  
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
                         ret <- scx %*% feature.loadings
                         return(ret)
                       }
                       ,mc.cores=ncores
    ) 
  cell.embeddings <- do.call(rbind, cell.embeddings)
  
  if(return_seurat_reduction){
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
                ,sd.pearson_residual = sd.pearson_residual))
      
  } else {
    return(list(cell.embeddings = cell.embeddings
                ,feature.loadings = feature.loadings
                ,mean.pearson_residual = mean.pearson_residual
                ,sd.pearson_residual = sd.pearson_residual))
  }
}
