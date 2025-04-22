

#' Bin cells into "plaids"
#'
#' @description R implementation of algorithm found in 'geosketch' python package
#' @param X feature matrix with cell IDs as rows and feature IDs as columns ()
#' @param N desired sample size
#' @param alpha divide feature space into uniform 'plaid' bins, where number of bins is in the range of
#' (`N * (1-alpha),  N * (1+alpha)`), where
#' @param max_iter maximum number of iterations used to achieve an acceptable
#'   minimum number of bins
#' @param return_grid determines whether or not to pass back bin labels instead of just the sample indices.
#'
geo_sketch_get_plaid <- function(X
                                 ,N = sqrt(nrow(X))
                                 ,alpha=0.05
                                 ,max_iter=200
                                 ,return_grid = FALSE
) {
  
  # Determine the total number of cells and compare it to the desired sample size 
  n_samples <- nrow(X)
  n_features <- ncol(X) 
  k <- N
  if (N > n_samples) {
    # Stop the function and return and error is the desired sample size is greater than the total number of cells given
    stop(paste0("N = ", N, " is greater than the number of cells in the feature matrix, n_samples = ", n_samples))
  }
  
  # Normalize features on a range between 0 and 1 in order to facilitate even binning across expression space
  X <- apply(X, 2, function(Y) (Y - min(Y)))
  X <- X/max(X)
  
  low_unit <- 0; high_unit <- 1
  unit <- (low_unit + high_unit) / 4
  
  # Iterate the number of bins each feature is broken into until the total
  # number of bins containing cells is larger than `(1-alpha)*N`
  iter <- 1 # Set starting iteration
  while (iter <= max_iter) { # Break loop when the number of iterations surpasses the defined maximum
    #message(paste0("Iteration Number ", as.character(iter))) # Report the current iteration the loop is on
   
    bins <- unique(c(seq(0, 1, by=unit),1)) # Define the bin ranges based on the current iteration
    #bins <- seq(0, 1, length.out=nbins) # Define the bin ranges based on the current iteration
    message(paste0("unit=", unit, ", bins: ", paste0(bins, collapse="--")))
    #message(paste0("unit=", 1/(nbins-1), ", bins: ", paste0(bins, collapse="--")))
    Xbins <- apply(X, 2, function(Y) .bincode(Y, bins, TRUE, TRUE)) # Bin cells across each feature using the given bin ranges
    
    griddt <- data.table(Xbins)[,grp:=.GRP,by=c(colnames(X))]
    len_grid <- griddt[,max(grp)]
    message("number of grid cells: ", len_grid)
    if(len_grid > k * (1 + alpha)){
      ### too many grid cells, increase unit 
      low_unit <- unit
      unit <- (unit + high_unit) / 2
      #nbins <- nbins + 1
    } else if (len_grid < k * (1-alpha)) {
      ### not enough grid cells, decrease unit 
      high_unit <- unit
      unit <- (unit + low_unit) / 2
    } else {
      break
    }
    iter <- iter + 1
  } 

  griddt[,idx:=1:.N] 
  griddt[,ncells_grp:=.N,by=grp] 
   
  cells_per_plaid <- griddt[,.N,by=.(grp)][,summary(N)]
  message("Summary of cells per plaid bin")
  print(cells_per_plaid)
  
  rand_samp <- griddt[sample(1:.N, N)]
  message("Summary of cells bin from uniform random sample")
  print(rand_samp[,.N,by=.(grp)][,summary(N)])
  message(rand_samp[,paste0("unique bins covered from random sample: ", uniqueN(grp))]) 
   
  samp_idx <- griddt[sample(1:.N)][,head(.SD,1),by=grp]
 
   
  if(nrow(samp_idx) < N){
    N_need <- N - nrow(samp_idx)
    samp_idx <- rbindlist(list(
          samp_idx
          ,griddt[!idx %in% samp_idx][sample(1:.N)][,head(.SD,1),by=.(grp)][1:N_need]
      ))
  }
  if(nrow(samp_idx) > N){
    #### Drop the least populated bins
    samp_idx <- samp_idx[order(-ncells_grp)][1:N]
  }
  return_idx <- sort(samp_idx[["idx"]])
  if(return_grid){
    return(list(griddt=griddt, idx = return_idx))
  } else {
    return(list(idx = return_idx))
  }
}
