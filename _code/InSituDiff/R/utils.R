
#' check format of x and obj
#' @param x Data matrix input to other functions
#' @param obj Output of initializeISD
#' @return Errors if needed
checkxandobj <- function(x, obj) {
  if (is.null(rownames(x))) {
    stop("x needs rownames")
  }
  if (is.null(colnames(x))) {
    stop("x needs colnames")
  }
  if (!identical(rownames(x), rownames(obj[[1]]))) {
    stop("rownames of x and the initializeISD obj don't match")
  }
}


#' Format "cells" argument to take the form of names
#'
#' @param cells Cells argument input to a larger function
#' @param x Data matrix input into a larger function
#' @return Cells formatted as a vector of cell IDs
formatCellsArg <- function(cells, x) {
  if (!is.null(cells)) {
    # if indices are input, convert to rownames:
    if (is.numeric(cells)) {
      cells <- rownames(x)[cells]
    }
    # if a logical vector is input, convert to indices
    if (is.logical(cells)) {
      cells <- which(cells)
    }
  } else {
    cells <- rownames(x)
  }
  return(cells)
}


#' Format "genes" argument to take the form of names
#'
#' @param genes genes argument input to a larger function
#' @param x Data matrix input into a larger function
#' @return Cells formatted as a vector of cell IDs
formatGenesArg <- function(genes, x) {
  if (!is.null(genes)) {
    # if indices are input, convert to rownames:
    if (is.numeric(genes)) {
      genes <- rownames(x)[genes]
    }
    # if a logical vector is input, convert to indices
    if (is.logical(genes)) {
      genes <- which(genes)
    }
    # if a list, assume it's modules, and format correctly:
    if (is.list(genes)) {
      # make sure each module gets a name:
      if (is.null(names(genes))) {
        names(genes) <- sapply(genes, function(x) {
          paste0(x[1:min(3,length(x))], collapse = "_")
        })
      }
    }
  } else {
    genes <- colnames(x)
  }
  return(genes)
}


#' Get a random-like but deterministic subset of a vector
#' 
#' @param vec Vector to subset
#' @param n Number of elements
#' @return A subset of the vector
#' @importFrom digest digest
pseudoRandomSample <- function(vec, n) {
  # hashes <- sapply(vec, function(x) as.numeric(paste0("0x", substr(digest::digest(x, algo = "xxhash64"), 1, 8))))
  # return(vec[order(hashes)[seq_len(min(n, length(hashes)))]])
  samp <- withr::with_seed(seed = 0, {use <- sample(seq_len(length(vec)), min(n, length(vec)), replace = FALSE)})
  return(samp)
}



