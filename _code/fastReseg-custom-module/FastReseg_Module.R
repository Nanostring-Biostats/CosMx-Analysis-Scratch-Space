# Name: FastReseg RNA Custom Module
message("FastReseg Custom Script Version: 1.1.2")

# Copyright 2023-2025 Bruker Spatial Biology, Inc.
# This software and any associated files are distributed pursuant to the Bruker Spatial Biology AtoMx Spatial
# Informatics Platform Software as a Service Agreement, available on the Bruker Spatial Biology, Inc
# website at www.nanostring.com, as updated.  All rights reserved.  No permission is granted to modify,
# publish, distribute, sublicense or sell copies of this software.

# Description: Evaluate cell segmentation error and trim off contaminating transcripts based on the spatial profile 
#              of RNA transcripts using FastReseg package. This module would first check the existence of FastReseg 
#              outcomes and run the evaluation if no existing data or different input hash. The FastReseg evaluation 
#              outcomes are saved as new soma called `RNA_trimmed` in current study. User can specify whether to 
#              overwrite the `RNA` soma with the new soma for reruning AtoMx pipeline from the step of QC  module. 
#              This module also adds a column called `lrtest_nlog10P` to the `RNA` obs and can be used in DE analysis. 

# User Defined Variables
# Configurations on setting up FastReseg evaluation:
# cellType     		  - cell type method, default = "InSituType_clusters" for InSituType, use "Leiden_clusters" for  
#                     leiden clustering or "NULL" to rely on external reference only; If cellType is not NULL and use_refProfiles 
#                     is selected, the shared clusters between the two will be used. 
# refProfiles_file 	- (optional) file path to external reference profiles (default to leave as blank); use Custom Module UI 
#                     to upload the reference matrix first if desired; if leave blank, use `cellType` only; If cellType is 
#                     not NULL and refProfiles_file has been uploaded, the shared clusters between the two will be used.  
# transNum_cutoff 	- minimal number of transcript per cell for segmentation evaluation (default = 50); in range of [3, Inf), 
#                     use 20 for 100plx, 50 for 1K-plx, 6K-plx or WTx panels. 
# flagCell_cutoff	  - minimal degree of spatial dependency for cells with contamination (default = 5); in range of [2, 10], 
#                     lower cutoff would flag more cells with putative contamination. 
# badScore_cutoff	  - maximum score for transcript with bad fit in current cell type (default = 2); in range of [1, 3], 
#                     lower cutoff is more permissive in flagging transcript-level contamination. 
# svm_gamma 	      - gamma used in svm function to separate good- vs. bad-fit transcripts (default = 0.4); in range of (0, 1], 
#                     smaller value for more spatial localized impact of bad-fit transcripts and thus more effective isolation 
#                     of bad-fit transcript groups with less smooth boundary.
# svm_cost	        - cost used in svm function to separate good- vs. bad-fit transcripts (default = 1); in range of (0, Inf], 
#                     higher value for more smooth spatial boundary with lower relative impact from singlets bad-fit transcripts.  
# maxCoreNum 		    - maximum core number allowed for parallel computing (default = 8); in range of [1, 96], use lower core 
#                     number for study with high transcript number per FOV. 

# Options on export or integrate FastReseg outcomes:
# overwrite_RNA_soma  - boolean to overwrite existing RNA soma with FastReseg outcomes (default = FALSE); if TRUE, all existing 
#                       post-QC data would be lost; recommend to export current data to seurat before overwriting. 

# Load Packages
# nanopipeline, FastReseg

# Module Code
library(nanopipeline)
library(FastReseg)
library(ggplot2)

# check if the custom module is built correctly ----
# Test User Variables
variableTest <- function(varName, varType, msg, required = TRUE){
  if(!varType %in% c("character", "logical", "numeric", "file")){
    stop("varType must be \"character\", \"logical\", \"numeric\", or \"file\"")
  }
  
  varMsg <- NULL
  
  typeClass <- switch(varType, 
                      character={"STRING"},
                      logical={"BOOL"},
                      numeric={"NUM"},
                      file={"FILE"})
  
  if(required == TRUE){
    if(!exists(varName)){
      varMsg <- paste0("\n\"", varName, "\" ", "was not set in custom module creation")
    }
  }
  
  if(exists(varName)){
    if(varType == "logical"){
      if(is.null(get(varName))){
        assign(varName, FALSE, envir = .GlobalEnv)
      }
    }
    
    if(varType == "file"){
      if(!is(get(varName), "character")){
        if(!is.null(get(varName)) | required == TRUE){
          varMsg <- paste0(varMsg, paste0("\n\"", varName, "\" varType was not set as ", typeClass, 
                                          " in custom module creation"))
        }
      }
      if(!file.exists(get(varName))){
        varMsg <- paste0(varMsg, paste0("\n\"", varName, "\" was not uploaded properly, please re-upload"))
      }
    }else{
      if(!is(get(varName), varType)){
        if(!is.null(get(varName)) | required == TRUE){
          varMsg <- paste0(varMsg, paste0("\n\"", varName, "\" varType was not set as ", typeClass, 
                                          " in custom module creation"))
        }
      }
    }
    
    
    if(required == TRUE){
      if(is.null(get(varName))){
        varMsg <- paste0(varMsg, paste0("\n\"", varName, "\" was not set and is required"))
      }else{
        if(get(varName) == ""){
          varMsg <- paste0(varMsg, paste0("\n\"", varName, "\" was not set and is required")) 
        }
      }
    }
  }
  
  if(!is.null(varMsg)){
    varMsg <- paste0(msg, varMsg)
  }else{
    varMsg <- msg
  }
  
  return(varMsg)
}


variableMsg <- NULL
variableMsg <- variableTest(varName = "cellType", varType = "character", required = TRUE, msg = variableMsg)
variableMsg <- variableTest(varName = "transNum_cutoff", varType = "numeric", required = TRUE, msg = variableMsg)
variableMsg <- variableTest(varName = "flagCell_cutoff", varType = "numeric", required = TRUE, msg = variableMsg)
variableMsg <- variableTest(varName = "badScore_cutoff", varType = "numeric", required = TRUE, msg = variableMsg)
variableMsg <- variableTest(varName = "svm_gamma", varType = "numeric", required = TRUE, msg = variableMsg)
variableMsg <- variableTest(varName = "svm_cost", varType = "numeric", required = TRUE, msg = variableMsg)
variableMsg <- variableTest(varName = "maxCoreNum", varType = "numeric", required = TRUE, msg = variableMsg)

variableMsg <- variableTest(varName = "overwrite_RNA_soma", varType = "logical", required = TRUE, msg = variableMsg)

# stop on variable errors if applicable
if(!is.null(variableMsg)){
  stop(variableMsg)
}

# check optional file, set use_refProfiles to TRUE if file exists 
use_refProfiles <- FALSE
refFileMsg <- NULL
refFileMsg <- variableTest(varName = "refProfiles_file", varType = "file", required = TRUE, msg = refFileMsg)
if(exists("refProfiles_file")){
  # file path provided
  if(!is.null(refFileMsg)){
    message("Current study has no external reference profiles uploaded via FastReseg module. Please upload reference `profile_matrix` (csv, rda or RData file) first if want to use it in segmentation evaluation. ")
    message("If don't want to use external reference profiles, please leave the file upload as blank.")
    stop(refFileMsg)
  }else{
    use_refProfiles <- TRUE # The file has uploaded. 
  }
}


if(flagCell_cutoff>10 || flagCell_cutoff <2 || badScore_cutoff <1 || badScore_cutoff >3 || transNum_cutoff <3){
  stop(sprintf("`flagCell_cutoff` = %d, `badScore_cutoff` = %d, and `transNum_cutoff` = %d, outside the allowed range of [2,10], [1,3] and [3, Inf) respectively. ", 
               flagCell_cutoff, badScore_cutoff, transNum_cutoff))
}

if(svm_gamma>1 || svm_gamma <=0 || svm_cost <=0){
  stop(sprintf("`svm_gamma` = %.3f, and `svm_cost` = %.3f, outside the allowed range of (0,1] and (0, Inf) respectively. ", 
               svm_gamma, svm_cost))
}

if(maxCoreNum <1){
  stop(sprintf("`maxCoreNum` = %d is not allowed. ", maxCoreNum))
}else{
  maxCoreNum <- ceiling(maxCoreNum)
}



#### functions ----
## required nanopipeline v1.2 functions:
# usableCores(), readInTranscriptCoords(), writeTranscriptCoords()

## nanopipeline v1.3 fastReseg functions:
# evalFun_perFOV_fastReseg(), runTileDBEvalFastReseg()
# internal method for tiledbsc::AnnotationDataframe: from_dataframe_batched()


## supporting functions: 
# check alignment and optional to extract transcript coords only
checkAndGetTranscriptCoords <- function(tiledbsc_dataset, 
                                        soma_names = c("negprobes", "falsecode"), 
                                        extract = TRUE){
  if(!is(tiledbsc_dataset, "SOMACollection")){
    stop("tiledbsc_dataset must be a SOMACollection")
  }
  
  soma_names <- intersect(soma_names, names(tiledbsc_dataset$somas))
  if(length(soma_names)<1){
    stop('Cannot find any of the input soma names in tildebsc_dataset.')
  }
  
  # initialize, return data.frame if extract = TRUE, return logic vector if extract = FALSE
  if(extract){
    combined_transCoords <- NULL
  } else {
    combined_transCoords <- c()
  }
  
  for(tgrt_soma in soma_names){
    if("transcriptCoords" %in% names(tiledbsc_dataset$somas[[tgrt_soma]]$obsm$members)){
      transCoord_colns <- tiledbsc_dataset$somas[[tgrt_soma]]$obsm$members[["transcriptCoords"]]$attrnames()
      if(any(!c('x_slide_mm','y_slide_mm') %in% transCoord_colns)){
        stop(paste0('Transcript coordinates have not been aligned for ', tgrt_soma, 
                    '; please create a new study with AtoMx v1.3 or higher.'))
        
      }
      if(extract){
        combined_transCoords <- rbind(combined_transCoords, 
                                      readInTranscriptCoords(tiledbsc_dataset$somas[[tgrt_soma]]$obsm$members$transcriptCoords))
      } else {
        combined_transCoords <- c(combined_transCoords, setNames(TRUE, nm = tgrt_soma))
      }
    }
  }
  
  return(combined_transCoords)
}


## supporting functions to extract config values:
# txt file to config value list
txtToConfig <- function(file_path){
  f <- readLines(file_path)
  # Create an empty list
  configVal_list <- list()
  
  # Iterate over the lines and extract the elements
  # not working for nested list or data.frame 
  for (line in f) {
    line <- gsub('\"', '', line, fixed = TRUE)
    
    if(line != ""){
      if(grepl("^\\$", line)){
        current_key <- gsub("^\\$", "", line)
        configVal_list[[current_key]] <- c()
        next
      }else if(grepl("^\\[\\d+\\] ", line)) {
        current_val <- gsub("^\\[\\d+\\] ", "", line)
        if(current_val != ""){
          current_val <- strsplit(current_val, " ")[[1]]
          current_val <- current_val[current_val!=""]
          # check if numeric 
          suppressWarnings( tmpNum <- as.numeric(current_val))
          if(!any(is.na(tmpNum))){
            current_val <- tmpNum
          } else {
            # check if logic 
            suppressWarnings(tmpFlag <- as.logical(current_val))
            if(!any(is.na(tmpFlag))){
              current_val <- tmpFlag
            }
          }
        }
        
      }else{
        # value directly, like NULL, character(0), or matrix (truncated)
        # return as it is after parse, won't handle matrix 
        tryCatch({
          current_val <- eval(parse(text = line))
        }, error = function(cond){
          current_val <- line
        })
        
      }

      configVal_list[[current_key]] <- c(configVal_list[[current_key]], current_val)
    }
    
  }
  return(configVal_list)
}

# aws config to local folder and then extract to value list 
awsToConfig <- function(S3Path, arry){
  s3_get(S3Path = S3Path, arry = arry, localPath = "/tmp/tmp.txt", fileType = "file")
  configVal_list <- txtToConfig(file_path = "/tmp/tmp.txt")
  file.remove("/tmp/tmp.txt")
  return(configVal_list)
}


# get file paths for config files, new function does NOT allow order by time, latest to older 
extractPathToConfigs <- function(studyDirectory){
  ## first assume s3 tileDB format with `//configs/` folder 
  configFiles <- tryCatch({
    output <- s3_ls(paste0(studyDirectory, "/configs/"), arry = study, fileType = "folder")
    output  # Return the output if successful
  }, error = function(e) {
    NA  # Return NA if an error occurs
  })
  if(all(is.na(configFiles))){
    configFiles <- data.frame(fullPath = character(0), 
                              fileName= character(0))
  }else{
    configFiles <- data.frame(fullPath = configFiles, 
                              fileName = gsub(paste0(studyDirectory, "/configs/"), "", configFiles))
    rownames(configFiles) <- gsub("\\.txt", "", configFiles$fileName)
  }
  
  ## further check 2nd place `/configs/` folder for config_loading if not found
  if(is.na(configFiles['config_loading', 'fullPath'])){
    cf2 <- tryCatch({
      output <- s3_ls(paste0(studyDirectory, "configs/"), arry = study, fileType = "folder")
      output  # Return the output if successful
    }, error = function(e) {
      NA  # Return NA if an error occurs
    })
    if(!all(is.na(cf2))){
      cf2 <- s3_ls(paste0(studyDirectory, "configs/"), arry = study, fileType = "folder")
      cf2 <- data.frame(fullPath = cf2, 
                        fileName = gsub(paste0(studyDirectory, "configs/"), "", cf2))
      rownames(cf2) <- gsub("\\.txt", "", cf2$fileName)
      
      cfNames <- setdiff(rownames(cf2), rownames(configFiles))
      if(length(cfNames)>0){
        configFiles <- rbind(configFiles, cf2[cfNames, , drop = FALSE])
      }
      rm(cfNames)
    }
    
    rm(cf2)
  }
  
  return(configFiles)
}

# add method  
tiledbsc::AnnotationDataframe$set("public", "from_dataframe_batched", overwrite = TRUE, 
                                  function(x, index_col = "obs_id", chunk_size = 1e6){
                                         nanopipeline::from_dataframe_batched(self, annots=x, index_col, chunk_size)
                                  })

study <- tiledbsc::SOMACollection$new(uri = studyDirectory, ctx = study$ctx, verbose = FALSE)

# check input arguments and initialize ----
if(!endsWith(studyDirectory, "/")){
  studyDirectory<- paste0(studyDirectory, "/")
}

configFiles <- extractPathToConfigs(studyDirectory)

# get config_loading 
if(is.na(configFiles['config_loading', 'fullPath'])){
  stop("Cannot find `config_loading` in the `configs` folder of current study: \n", studyDirectory)
}
config_loading <- awsToConfig(S3Path = configFiles['config_loading', 'fullPath'], 
                              arry = study)

# check cell Type 
cellType <- match.arg(cellType, c("NULL", "InSituType_clusters", "Leiden_clusters"))

if(cellType == "NULL"){
  cellType <- NULL
  message("CellType = NULL, use the external reference profiles for segmenation evaluation.")
  if(!use_refProfiles){
    stop(refFileMsg) # refProfiles_file not uploaded
  }
  
} else {
  obs_colns <- study$somas$RNA$obs$attrnames()
  if(cellType == "InSituType_clusters"){
    # guess by the column names
    # column on the left is the latest for insitutype
    cellType <- grep("^RNA_.*_posterior_probability$", obs_colns, value = TRUE)[1]
    cellType <- gsub("_posterior_probability$", "_clusters", cellType)

    ## check config file if not found 
    if("config_nbclust" %in% rownames(configFiles) & is.na(cellType)){
      config_nbclust <- awsToConfig(S3Path = configFiles['config_nbclust', 'fullPath'],
                                    arry = study)
      cellType <- paste(config_nbclust$assay,config_nbclust$seurat_slot,"clusters",sep="_")
      rm(config_nbclust)
    }
    

  }else{
    # column on the right is the latest for leiden
    cellType <- utils::tail(grep("^nn_.*_cluster_cluster_", obs_colns, value = TRUE), 1)
    
    ## check config file if not found 
    if("config_cluster" %in% rownames(configFiles) & length(cellType)<1){
      config_cluster <- awsToConfig(S3Path = configFiles['config_cluster', 'fullPath'],
                                    arry = study)
      cellType <- paste(config_cluster$seurat_slot, "cluster", config_cluster$cluster_name, sep="_")
      rm(config_cluster)
    }
  }
  
  cellType <- intersect(cellType, obs_colns)
  
  if(length(cellType)<1){
    stop("The provided cellType does not exist in current data, check if Cell Typing (InSituType) or Leiden Clustering module has been ran.")
  }
}

# reference 
reference_profiles <- NULL
if(use_refProfiles){
  # read profiles from uploaded file 
  message("Load external reference from uploaded file")
  
  # read in and pass to config, won't be recoverable from txt but only compare hash here
  if(grepl('\\.csv$|\\.txt$', refProfiles_file, ignore.case = TRUE)){
    reference_profiles <- read.csv(refProfiles_file, row.names = 1, header = TRUE, 
                                   sep = ifelse(grepl("\\.csv$", refProfiles_file, ignore.case = TRUE), ",", "\t"))
    reference_profiles <- as.matrix(reference_profiles)
    
  }else if(grepl('\\.RData$|\\.rda$', refProfiles_file, ignore.case = TRUE)){
    load(refProfiles_file)
    if(!exists("profile_matrix")){
      stop("External reference profiles must be stored under `profile_matrix` variable in uploaded .RData or .rda file.")
    }
    reference_profiles <- as.matrix(profile_matrix)
    rm(profile_matrix)
  }else{
    stop("reference_profiles must be a .csv, .txt or .RData, .rda file")
  }
} 

if(is.null(cellType) && is.null(reference_profiles)){
  stop("Must provide either `cellType` or external reference profiles. Please use Cell Type (InSituType) fundation module to upload external reference matrix first if use_refProfiles = TRUE.")
}

config_evalFastReseg <- list(
  seurat_slot = "RNA_trimmed", # output soma name
  assay = "RNA", 
  cellType = cellType, 
  reference_profiles = reference_profiles, # actual ref matrix
  flagModel_TransNum_cutoff = transNum_cutoff, 
  flagCell_lrtest_cutoff = flagCell_cutoff, 
  svmClass_score_cutoff = 0 - badScore_cutoff,  
  svm_scale = FALSE, 
  svm_gamma = svm_gamma, 
  svm_cost = svm_cost, 
  remove_new_data =  TRUE # always true for tiledbsc due to requirement of same cells 
)

## check if FastReseg data exist ----
runEval_flag <- TRUE
overwrite <- FALSE
if(config_evalFastReseg$seurat_slot %in% names(study$somas)){
  oldConfigHash <- study$somas[[config_evalFastReseg$assay]]$get_metadata(prefix = "evalFastReseg_config_parameters")
  if(verifySeuratConfigHash(config_evalFastReseg, oldConfigHash)){
    print("Found exisiting `RNA_trimmed` soma and matched configurations, skip the evaluation module. ")
    runEval_flag <- FALSE 
  } else {
    print("Existing `RNA_trimmed` soma was found but the configurations are not matched with input.")
    print("Change overwrite to TRUE, erase existing `RNA_trimmed` soma to redo evaluation step...")
    overwrite <- TRUE
    
  }
}


## run evaluation module ----
if(runEval_flag){
  # check alignment of RNA transcript coords, error if it hasn't been done. 
  alignRes <- checkAndGetTranscriptCoords(tiledbsc_dataset = study, 
                                          soma_names = "RNA", 
                                          extract = FALSE)
  
  nanopipeline:::write_config(scdataset = study, outFolder = studyDirectory, 
                              config = config_evalFastReseg, config_name = "config_evalFastReseg",
							  s3_uri = studyDirectory)
  
  # check bigData for earlier abortion 
  bigData <- nanopipeline:::checkBigData(study$somas[[config_evalFastReseg$assay]], threshold = 2^28)
  if(bigData){
    stop("Current FastReseg evaluation does not support big dataset.")
  }
  
  # cap the core number 
  options(mc.cores = maxCoreNum)  
  
  runTileDBEvalFastReseg(config_evalFastReseg = config_evalFastReseg,
                         config_loading = config_loading,
                         scdataset = study,
                         tiledb_path = studyDirectory,
                         overwrite = overwrite)
  options(mc.cores = NULL)
    
  # Update new additions to somas for tiledb_scdataset to maintain
  study$initialize(studyDirectory)
}

rm(reference_profiles)

## create summary plot for resegmentation outcomes ----
outDir <- "/output"
dir.create(outDir)
obs <- study$somas[['RNA']]$members$obs$to_dataframe()
if(!all(c("lrtest_nlog10P", "flagged") %in% colnames(obs))){
  stop("No FastReseg evaluation has been run on current study.")
}
fig <- ggplot2::ggplot(obs, ggplot2::aes(x = lrtest_nlog10P, fill = as.factor(flagged)))+
  ggplot2::geom_histogram()+
  ggplot2::labs(title = sprintf("Total %d cells: %d cells flagged by FastReseg.",  
                                nrow(obs),  sum(obs$flagged, na.rm = TRUE)),
                subtitle = sprintf("%d cells below evaluation cutoff.",  
                                sum(is.na(obs$flagged))), 
                fill = 'flagged')


tryCatch(
  {
    ggplot2::ggsave(
      filename = fs::path(
        outDir,
        "FastReseg_segmentation_flagging_score_hist.png"
      ),
      plot = fig,
      width = 6,
      height = 6,
      units = "in",
      dpi = 150
    )
  },
  error = function(e) {
    message("ggsave PNG failed: ", conditionMessage(e))
  }
)

rm(obs, fig)
gc()



## overwrite RNA soma with RNA_trimed ----
if(overwrite_RNA_soma){
  print("Overwrite current RNA soma with new RNA_trimmed soma. All post-QC data like dimension reduction, DE results would be lost.")
  print("Once RNA soma is overwritten, one should rerun the pipeline from QC module.")
  
  # write latest.fovs to RNA_trimmed
  if("latest.fovs" %in% names(study$somas[["RNA"]]$obsm$members)){
    if(! "latest.fovs" %in% names(study$somas[[config_evalFastReseg$seurat_slot]]$obsm$members)){
      latest.fovs <- study$somas[["RNA"]]$obsm$members$latest.fovs$to_matrix()
      study$somas[[config_evalFastReseg$seurat_slot]]$obsm$add_annotation_matrix(latest.fovs, 
                                                                                 "latest.fovs")
      rm(latest.fovs)
    }
  }
  
  # check and update transcript coords for negprobes and falsecodes
  alignRes <- checkAndGetTranscriptCoords(
    tiledbsc_dataset = study, 
    soma_names = c("negprobes", "falsecode"),  
    extract = FALSE)
  
  # add in new columns for soma with transcriptCoords
  for(tgrt_soma in names(which(alignRes))){
    tCoordsDF <- readInTranscriptCoords(study$somas[[tgrt_soma]]$obsm$members[['transcriptCoords']])
    if(!"z_slide_mm" %in% colnames(tCoordsDF)){
      tCoordsDF[["z_slide_mm"]] <-  tCoordsDF[["z_FOV_slice"]]*0.8*0.001
    }
    tCoordsDF[["trimmed"]] <- 0L
    tCoordsDF[["old_cell_id"]] <- tCoordsDF[["cell_id"]]
    writeTranscriptCoords(study$somas[[tgrt_soma]]$obsm, tCoordsDF)
    rm(tCoordsDF)
  }
  
  
  rna_soma_uri <- study$somas$RNA$uri
  trimmed_soma_uri <- study$somas[[config_evalFastReseg$seurat_slot]]$uri
  
  # delete old RNA soma
  manageTileDBStoredData(study, rna_soma_uri, overwrite = TRUE)
  
  # initialize place holder
  soma_rna <- tiledbsc::SOMA$new(rna_soma_uri, verbose = FALSE, ctx = study$ctx)
  
  # move RNA_trimmed to RNA soma
  s3_sync(inPath = paste0(trimmed_soma_uri, '/'), 
          arry = study$somas[[config_evalFastReseg$seurat_slot]],
          outPath = paste0(rna_soma_uri, '/'), 
          fileType = "folder")
  s3_rm(S3Path = paste0(trimmed_soma_uri, '/'), 
        arry = study$somas[[config_evalFastReseg$seurat_slot]], 
        fileType = "folder", recursive = TRUE)

  study$remove_member(config_evalFastReseg$seurat_slot)
  
  # update study with soma change 
  study$initialize(studyDirectory)
  
  print("RNA soma now contains the post-trimmed data after FastReseg evaluation. All preivous post-QC data is lost except for the ones in `obs`.")
  
  rm(rna_soma_uri, trimmed_soma_uri, soma_rna)
}
