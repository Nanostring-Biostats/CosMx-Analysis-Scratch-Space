
#' Project PCA onto new data 
#' 
#' @description
#' 
#' Take the results generated from `sparse_quasipoisson_pca_seurat` and project them onto new cells.
#' 
#'
#' @param x A sparse genes x cells counts matrix 
#' @param pcaobj generated from `sparse_quasipoisson_pca_seurat`  
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
                                                   ,totalcounts = NULL
                                                   ,ncores=1
                                                   ,verbose=TRUE
                                                   ,block.size=10**4
                                                   ,ndims.print = 1:5
                                                   ,nfeatures.print = 30
                                                   ,return_seurat_reduction = TRUE
){
 
  x <- check_x_is_dgcmatrix(x)
  
  nn <- ncol(x)
  phi <- check_phi(x, pcaobj$phi) # quasi_poisson_variance_inflation   
  stopifnot(length(phi) %in% c(1, length(pcaobj$grate)))
  totalcounts <- check_totalcounts(x, totalcounts)
  grate <- pcaobj$grate
  stopifnot(!is.null(names(grate)))
  stopifnot(all(names(grate) %in% rownames(x))) 
  x <- x[names(grate),] 
  grate <- check_grate(x, grate, totalcounts)
  
  ### diagonal matrix (1/(phi*estimated gene frequency))
  root_grate_phi_diag <- Matrix::Diagonal(x = sqrt(1/(grate * phi)))
  dimnames(root_grate_phi_diag) <- list(names(grate) , names(grate))
  
  ### diagonal matrix (1/totalcounts)
  root_tc_diag <- Matrix::Diagonal(x = sqrt(1/totalcounts))
  dimnames(root_tc_diag) <- list(colnames(x) , colnames(x))
  
  ### y / sqrt(v(y)) 
  ytilde <- root_grate_phi_diag %*% x %*% root_tc_diag ## sparse  n x g
  
  ### gene-wise average of y/sqrt(v(y)) - muhat/sqrt(v(y))
  ### Note: muhat / sqrt(v(y)) = muhat/sqrt(muhat*phi) = sqrt(muhat/phi)
  mean.pearson_residual <- pcaobj$mean.pearson_residual
  
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
#    cellnames_max_val_vec <- colnames(ytilde)[findInterval(seq(ytilde@x)-1, ytilde@p[-1]) + 1]
    cellnames_max_val_vec <- colnames(ytilde)[rep.int(seq_len(ncol(ytilde)), diff(ytilde@p))]
    
    max_val_vec <- max_val_vec + 
      sqrt(grate[names(max_val_vec)]/phi) * sqrt(totalcounts[cellnames_max_val_vec])
    xvec <- pmin(xvec, max_val_vec) 
    
    ytilde@x <- xvec
    rm(list=c("xvec", "max_val_vec", "cellnames_max_val_vec")); gc()
  }
  
  #### Get the eigenvectors / loadings 
  feature.loadings <- pcaobj$reduction.data@feature.loadings
  feature.loadings <- feature.loadings[rownames(x),]
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
