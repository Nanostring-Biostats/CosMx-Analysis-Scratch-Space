
library(reticulate)
library(data.table)
library(Matrix)
library(Seurat)
library(plyr)
library(dplyr)
library(ggplot2)
library(ggExtra)
library(kableExtra)

# Simple function render with include_graphics. Copies file to assets folder if not there.
render <- function(target_parent_dir, filename, source_parent_folder=NULL, overwrite=FALSE, verbose=FALSE){
  
  if(!dir.exists(target_parent_dir)){
    dir.create(target_parent_dir, recursive = TRUE)
  }
  
  target_path <- file.path(target_parent_dir, filename)
  
  if(!is.null(source_parent_folder)){
    source_path = file.path(source_parent_folder, filename)
    if(!file.exists(source_path)){
      stop("source_parent_folder is not null but it also does not exists. check path name.")
    }
    if(!file.exists(target_path) | overwrite){
      file.copy(from = source_path, to = target_path, overwrite = TRUE)
    }
  }
  
  if(!file.exists(target_path)){
    stop("Could not find target file. check path name.")
  }
  return(knitr::include_graphics(target_path))
}


# simple function to put files in numeric order (not alphabetic order)
sortFiles <- function(x, i){
  cluster_number <- as.numeric(gsub(".png", "", unlist(lapply(strsplit(basename(x), split="_"), "[[", as.integer(i)))))
  out <- data.frame(file=x, cluster=cluster_number) %>% arrange(cluster) %>% pull(file)
  return(out)
}


