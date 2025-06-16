
#' Compute gene frequency from counts matrix. 
#' 
#' @description
#' 
#' Compute gene frequency (possibly by batch) , used in calculation of pearson residuals. 
#' 
#'
#' @param x A sparse genes x cells counts matrix 
#' @param obs A data frame of metadata information, containing columns for the 'batch variable' and 'cell id'.  
#' Or, row names of the data frame indicating 'cell id'.
#' @param batch_variable column in `obs` indicating 'batch'.  i.e., 'patient' or 'samplename' or 'slide_id'.
#' @param cellid_colname column in `obs` indicating the cell identifier, i.e., 'cell_ID'.
#'
#'
#' @export  
gene_frequency <- function(x
                           ,obs = NULL
                           ,batch_variable =NULL
                           ,cellid_colname = "cell_ID"
){
  
 
 
  if(is.null(obs) || is.null(batch_variable)){
    if(!is.null(batch_variable)){
      warning(paste0("ignoring batch_variable='",batch_variable, "', no obs dataframe of metadata specified"))
    }
    if(!is.null(obs)){
      warning(paste0("ignoring argument passed for 'obs', no 'batch_variable' specified."))
    }
    grate <- Matrix::rowSums(x)
    if(any(is.na(grate))){
      stop("NA values found in these genes: ", paste0(names(grate[is.na(grate)]), collapse=","))
    }
    if(any(is.na(grate))){
      stop("Infinite values found in these genes: ", paste0(names(grate[is.na(grate)]), collapse=","))
    }
    grate <- grate / sum(grate)
  } else {
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
    
    ### expression rates by batch (\hat{p})
    grate <- x %*% batch_mat %*% Matrix::Diagonal(x=1/Matrix::colSums(batch_mat), names = colnames(batch_mat))
    grate <- grate %*% Matrix::Diagonal(x = 1/Matrix::colSums(grate),names=colnames(batch_mat))
    
  } 
  return(grate) 
}
