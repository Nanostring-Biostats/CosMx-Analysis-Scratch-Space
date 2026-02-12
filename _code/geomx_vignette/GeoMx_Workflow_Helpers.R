#' @name classifyPAdjustMethod
#' @rdname classifyPAdjustMethod
#' @title Classify Adjusted P-value Method
#' @description Given p.adjust method, classify method to either: FDR, FWER or P.
#' @param \code{p_adj_method} A character specifying which multiple testing correction method to use. Options are those from \code{p.adjust.methods}. Defaults to "fdr"
#' @details The function intakes a \code{p.adjust.methods} character and outputs it's respective adjusted p-value class. The methods "fdr", "BH", and "BY" are classified as "FDR"; the methods "holm", "hochberg", "hommel", and "bonferroni" are classified as "FWER"; and "none" is classified as "P". Note that the input method can be of any case.
#' @export
#' @examples
#' classifyPAdjustMethod("FDR")
#' classifyPAdjustMethod(p_adj_method = "bonferroni")


classifyPAdjustMethod <- function(p_adj_method) {
    
 # ################
 # Check User Input
 # ################

 # Ensure the p_adj_method input is a character of length 1
 if(!inherits(p_adj_method, "character")){
  stop("Please check the input; input must be a string.")
 }
 if(!length(p_adj_method) == 1){
   stop("Please check the input; input must be a single string")
 }   

 # ########
 # Process
 # ######## 

    if(any(c("FDR", "BH", "BY") == toupper(p_adj_method))){
        p_adj_class <- "FDR"
    } else if (any(c("HOLM", "HOCHBERG", "HOMMEL", "BONFERRONI") == toupper(p_adj_method))) {
        p_adj_class <- "FWER"
    } else if (any(c("NONE") == toupper(p_adj_method))) {
        p_adj_class <- "P"
    } else {
        stop("Please select a valid p.adjust method")
    }

 # #######
 # Return
 # #######

 return(p_adj_class)
}

#' @name formatLMMResults
#' @rdname formatLMMResults
#' @title Format LMM Results from GeomxTools::mixedModelDE
#' @param lmm_results an object of class matrix that contains the LMM results from GeomxTools::mixedModelDE
#' @param p_adjust_method character. The method used to adjust p-values (see \code{p.adjust})
#' @details Note while the lsmeans is based on "contrast" (e.g., B - A; B minus A), the formatted data
#' that are returned are in terms of "Comparisons" (e.g., A vs B). In both cases, A refers to the base level.
#' @return a data.frame that is formatted with column names "Feature", "Comparison", "Estimate", "P", <toupper(p_adjust_method)>
#' @importFrom plyr ddply
#' @importFrom plyr .
#' @export
#' @examples
#' library(plyr)
#' library(GeomxTools)
#' data("kidney_targets")
#' kidney_targets_subset <- kidney_targets[1:5, ]
#' kidney_targets_subset$slide <- factor(kidney_targets_subset$slide)
#' kidney_targets_subset$pathology <- factor(kidney_targets_subset$pathology)
#' de <- mixedModelDE(kidney_targets_subset,
#'                   elt = "log_q",
#'                   modelFormula = ~ pathology + (1|slide),
#'                   groupVar = "pathology",
#'                   nCores = parallel::detectCores()-1,
#'                   multiCore = FALSE)
#' formatLMMResults(de)

formatLMMResults <- function(lmm_results, p_adjust_method = "fdr") {
 if(!inherits(lmm_results, "matrix")){
  stop("lmm_results needs be a matrix, See ?GeomxTools::mixedModelDE.")
 }
 if(!all(c("anova", "lsmeans") %in% rownames(lmm_results))){
  stop("Expected row names of lmm_results to have anova and lsmeans.")
 }
 df <- do.call(rbind, lmm_results["lsmeans", ])
 contrasts <- rownames(df)

 # Make sure there are not multiple " - " present in contrasts
 for(i in unique(contrasts)){
   if(length(strsplit(i, split=" - ")[[1]])>2){
     stop(paste0("Contrast \'", i, "\' has more than two split points. Please rename contrasts first."))
   }
 }
 # Contrast are of the for B - A. Convert to comparisons of the form
 # A vs B.
 df <- as.data.frame(df)
 contrast_pairs <- strsplit(contrasts, split=" - ")
 df$Comparison <- paste0(unlist(lapply(contrast_pairs, "[[", 2L)), " vs ", unlist(lapply(contrast_pairs, "[[", 1L)))
 colnames(df)[which(names(df) == "Pr(>|t|)")] <- "P"
 row.names(df) <- NULL

 # Add feature names
 df$Feature <- rep(colnames(lmm_results), each=nrow(lmm_results["lsmeans",][[1]]))

 # P-adjustment based on subsets of data faceted by Comparison
 df <- ddply(df, .(Comparison), function(x){
   x$padj <- p.adjust(x$P, method = p_adjust_method)
   return(x)
 })
 df <- df[, c("Feature", "Comparison", "Estimate","P", "padj")]
 colnames(df)[colnames(df)=="padj"] <- toupper(p_adjust_method) # 'FDR' used in standard Report


 return(df)
}

#' @name getColorPalette
#' @rdname getColorPalette
#' @title Set Color Palette
#' @description Parse and set colors to match annotations and be consistent throughout the report.
#' @param \code{input} a \code{data.frame} consisting of factors of interest with levels to be assigned colors
#' @param \code{start} what color in the palette to start with, default 1
#' @param \code{custom} Custom color palette to use instead of the predefined palette (a \code{vector} of mode \code{character})
#' @param \code{n} Number of colors needed
#' @param \code{method} Type of color palette to return: 'Main', 'Map', or 'Other'
#' @details  This function parses and sets colors to match annotations and be consistent throughout the report.
#' @return  A list of length equal to the number of columns in \code{input}. Names of the elements of the list are equal to the column names in \code{input}. Within each element of the list, a named character vector giving the discrete color value.
#'
#' @export
#' @examples
#' library(GeomxTools)
#' library(RColorBrewer)
#' color_pal <- getColorPalette(pData(sDAS::kidney_targets)[,c('class','segment','region')])
#' custom_pal <- c("#1B9E77", "#D95F02", "#7570B3", "#E7298A", "#66A61E", "#E6AB02", "#A6761D")
#' color_pal2 <- getColorPalette(pData(sDAS::kidney_targets)[,c('class','segment','region')], custom=custom_pal)
#' color_pal3 <- getColorPalette(pData(sDAS::kidney_targets)[,c('class','segment','region')], custom=brewer.pal(7,"Set1"))

getColorPalette <- function(input,              # Input must be data.frame
                            start = 1 ,         # Which color to start with
                            custom = NULL,      # A custom color palette, vector of type character
                            method = "Map"){    # Color palette type (this could change in the future)

## Check input
   # input should be a data.frame
   if(!inherits(input, "data.frame")){
      stop("input must be a data.frame")
   }

   # start must be a positive numeric integer.
   if(!is.null(start)) {
      if (!inherits(start, "numeric")) {
        stop("start must be numeric.")
      }
      if (start %% 1 != 0) {
        stop(paste0("The start given, ", start, " is not an integer."))
      }
      if (start < 1) {
        stop("start must be >= 1.")
      }
   }

   # custom must be a vector of mode character
   if(!is.null(custom)) {
      if(!inherits(custom, "character")){
        stop("custom color palette must be a vector of mode character")
      }
      if(length(custom) < 1){
        stop("custom must be at least length 1")
      }
      if(is.vector(custom) & is.list(custom)){
        stop("custom must be a vector of mode character, not a list")
      }
   }

   # method mustbe a character either: "Map", "Main", or "Other"; "Map" is default
   if(!method %in% c("Map","Main","Other")){
    stop("method must be one of: 'Main', 'Map', or 'Other'")
   }


# Make a palette data frame for reference at start of report
           l <- start
           input <- as.data.frame(input)
           color_pal <- list()
           input[] <- lapply(input, factor)
           for (i in 1:ncol(input)) {
            input[, i] <-
             factor(input[, i], levels = levels(factor(input[, i])), order = TRUE)
            pal <- report_pals(method = method,
                               n = length(levels(factor(input[, i]))),
                               start = l,
                               custom = custom)
            names(pal) <- levels(input[, i])
            l <- l + length(levels(input[, i]))
            color_pal[[length(color_pal) + 1]] <- pal
           }
           names(color_pal) <- names(input)

  return(color_pal)
}





#Helper function:
# Grab colors depending on which method: Main, Map, other
report_pals <- function(method = NULL,          # which palette type (Main, Map, or Other)
                        n = NULL,                   # how many colors to include in the palette
                        start = NULL,                # what color to start on
                        custom = NULL) {          # list of custom colors to be used as a palette
 if(!is.null(custom)) {
  pal <- custom
 } else {
  if(method == "Main") {
   # define palette for use in boxplots / scatter plots with annotations, n = 10 colors
   pal <- c("#3A6CA1", "#FFD861", "#CF4244", "#47BAB4", "#474747", "#EB739E", "#318026", "#A66293", "#F28E2B", "#8F6954")
   #defaults: blues  ,  yellows ,   reds   ,   teals  ,   grays  ,   pinks  ,   greens ,  purples ,  oranges ,  browns
  } else if(method == "Map") {
   # define palette to be passed to annotation maps (heatmap), n = 20 colors
   pal <- c("#3A6CA1", "#FFD861", "#D86769", "#AEE8E2", "#999999", "#FABFD2", "#318026", "#A66293", "#F28E2B", "#E3C0AC",
            "#A0CBE8", "#9E7E20", "#FFBFBD", "#2CABA3", "#474747", "#EB739E", "#A0E391", "#E8C1DE", "#FCB36A", "#B0846B")
  } else if(method == "Other") {
   # define palette to be passed to other types of factors that could be useful in the future: cell types, genes, etc, n = 20 colors
   pal <- c("#3A6CA1", "#FFD861", "#D86769", "#AEE8E2", "#999999", "#FABFD2", "#318026", "#A66293", "#F28E2B", "#E3C0AC",
            "#A0CBE8", "#9E7E20", "#FFBFBD", "#2CABA3", "#474747", "#EB739E", "#A0E391", "#E8C1DE", "#FCB36A", "#B0846B")
  }
 }

 ntot <- n + start - 1
 if(ntot > length(pal)) {
  pal <- rep(pal, ceiling(ntot / length(pal)))
 }
 return(pal[seq(start, ntot)])
}
#' @name getCV
#' @rdname getCV
#' @title Get the coefficient of variation
#' @description This function computes the coefficient of variation
#' (CV) of a vector of numeric or integer values. There must be at least
#' two non-NA values to calculate the CV.
#' @export
#' @examples
#' getCV(1:10)
#'

getCV <- function(vec){
  if(!inherits(vec, "numeric") & !inherits(vec, "integer")){
    stop("vec needs to be numeric.")
  }
  if(length(vec[!is.na(vec)])<2){
    stop("vec needs at least 2 non-NA observations.")
  }
  cv <- sd(vec)/mean(vec)
  return(cv)
}


setGeneric(name="getPCA", signature = c("object"),
           function(object, ...){
             standardGeneric("getPCA")
           })
#' @name getPCA
#' @rdname getPCA
#' @title Wrapper function for FactoMineR's PCA
#' @description Run FactoMineR's PCA with additional options.
#' @param \code{object} a NanoStringGeoMxSet object.
#' @param \code{elt} \code{character} specifying the expression \code{matrix}.
#' @param \code{use_top_cv} Logical. Use the targets with the high coefficient of variation? Default is FALSE.
#' @param \code{top_n} Numeric integer >=1. Can be NULL. If \code{use_top_cv} is TRUE, use this many targets. If \code{top_n}
#' greater than the number of features, no subsetting will be used.
#' @param \code{log2_transform} Logical. To log2 transform data or not. Default is FALSE.
#' @param \code{PCA_args} \code{pairlist} of arguments to pass to the PCA function. See details.
#' @param \code{the_prefix} \code{character} giving the prefix to use for column names. Default is "PCA".
#' @details This function is a wrapper for running FactoMineR's PCA function with a \code{NanoStringGeoMxSet}
#' object with custom and built-in (_i.e._, FactoMineR) arguments. For the built-in arguments,
#' \code{getPCA} passes the \code{pairlist} from the formal arguments of FactoMineR::PCA. Other than the 'X', 'graph', and 'axes' values,
#' the values in the \code{pairlist} will get passed to the PCA function. The user can optionally pass a \code{pairlist}
#'  object with different values. Note that the names of the element are required to be the same as the formal args version.
#'  Also note that no additional error or class checking is performed on those elements.
#'
#'  The PCA coordinates for individuals will be added as a new assayDataElement with the number of columns
#'  equal to the number specified in PCA_args$ncp (one row for each sample in object). The new element name will begin
#'  with \code{the_prefix}
#'
#' @return A modified \code{NanoStringGeoMxSet} object. The PCA coordinates for individuals will be added to pData(object)
#' with the number of new columns equal to the number specified in PCA_args$ncp.
#'
#' @importFrom Biobase pData
#' @importFrom Biobase fData
#' @importFrom Biobase esApply
#' @importFrom FactoMineR PCA
#' @import GeomxTools
#' @export
#' @examples
#' library(GeomxTools)
#' kidney_targets <- updateGeoMxSet(kidney_targets)
#' getPCA(kidney_targets, elt="log_q")


setMethod(f="getPCA",
          signature = "NanoStringGeoMxSet",
          definition = function(object,
                                elt,
                                use_top_cv=FALSE,
                                top_n=100,
                                log2_transform=FALSE,
                                PCA_args=formals(FactoMineR::PCA),
                                the_prefix="PCA") {

  ## Check the input

  # elt must be a non-zero length character and specifies an expression matrix
  if(!inherits(elt, "character")){
    stop("Please provide a character to elt.")
  }
  if(nchar(elt) == 0){
      stop("elt can not be a zero-length variable name.")
  }
  if(!inherits(assayDataElement(object, elt = elt), "matrix")){
    stop("elt was not found in object. Please check spelling.")
  }

  # use_top_cv must be logical.
  if(!is.logical(use_top_cv)){
    stop("use_top_cv must be logical.")
  }

  # top_n must be NULL or a positive numeric integer. If NULL, use_top_cv must be FALSE.
  if(is.null(top_n) & use_top_cv){
    stop("Please provide the number of targets to use when restricting to the top CV targets.")
  }
  if(!is.null(top_n)){
    if(!inherits(top_n, "numeric")){
      stop("top_n must be numeric.")
    }
    if(top_n %% 1 != 0){
      stop(paste0("The top_n given, ", top_n, " is not an integer."))
    }
    if(top_n <= 0){
      stop("top_n must be >= 1.")
    }
  }

  # log2_transform must be logical. If it is TRUE, no elements
  # in the expression matrix may be zero.
  if(!is.logical(log2_transform)){
    stop("log2_transform must be logical.")
  }
  if(log2_transform & any(assayDataElement(object, elt=elt)<=0)){
    stop("Cannot do log2 transformation on values less than or equal to zero.")
  }

  # PCA_args must be a pairlist with named elements: X, scale.unit, ncp
  # ind.sup, quanti.sup, quali.sup, row.w, col.w, graph, axes.
  if(!inherits(PCA_args, "pairlist")){
    stop("PCA_args must be a pairlist. See formals(FactoMineR::PCA).")
  }
  valid_names <- names(formals(FactoMineR::PCA))
  if(!all(valid_names %in% names(PCA_args))){
    missing_names <- valid_names[!valid_names %in% names(PCA_args)]
    stop(paste0("Elements in PCA_args must include all elements from formals(FactoMineR::PCA): ",
                paste(missing_names, collapse=", "), " were not found."))
  }

  # the_prefix must be a character string
  if(!inherits(the_prefix, "character")){
    stop("the_prefix must be a character.")
  }

  # New column names added (with the form <the_prefix>_coords_1, ..., <the_prefix>_coords_<PCA$ncp>)
  # must not already be present in the phenotype data.
  new_coord_names <- paste0(the_prefix, "_coords_",
                            as.character(seq.int(from=1, to=PCA_args$ncp, by=1)))
  if(any(new_coord_names %in% colnames(pData(object)))){
    stop("One or more coordinate names already exists. Please specify a different the_prefix value.")
  }

  # Analysis
  if(log2_transform){
    mat <- t(log2(assayDataElement(object, elt = elt)))
  } else {
    mat <- t(assayDataElement(object, elt = elt))
  }
  if(use_top_cv){
    mat <- getTopCV(mat, top_n)
  }
  pca <- FactoMineR::PCA(mat,
                         scale.unit = PCA_args$scale.unit,
                         ncp = PCA_args$ncp,
                         ind.sup = PCA_args$ind.sup,
                         quanti.sup = PCA_args$quanti.sup,
                         quali.sup = PCA_args$quali.sup,
                         row.w = PCA_args$row.w,
                         col.w = PCA_args$col.w,
                         graph=FALSE
                         )

  # Adjust column names.
  colnames(pca$ind$coord) <- gsub("Dim.",
                                  paste0(the_prefix, "_coords_"),
                                  colnames(pca$ind$coord))
  # Add columns of pca$ind$coord to phenotype data.
  pData(object) <-  cbind(
      as(pData(object), "data.frame"),
      as.data.frame(pca$ind$coord))

  # Append variation explained to experiment data
  pca_var_explained <- list(pca$eig[seq.int(from=1, to=PCA_args$ncp, by=1),"percentage of variance"])
  names(pca_var_explained) <- the_prefix
  experimentData(object)@other <- c(experimentData(object)@other, pca_var_explained)

  # return object
  return(object)

})

#' @title Helper function for selecting features with the greatest coefficient of variation
#' @param \code{mat} a matrix with columns as features and samples as rows
#' @param \code{top_n} an integer giving the number of features. If top_n >= ncol(mat), no subsetting is done.
#' @return a matrix with top_n columns.

getTopCV <- function(mat, top_n){
    if(top_n < ncol(mat)){
      cv <- apply(mat, 2, function(x){
        return(sd(x)/mean(x))
      })
      top_cv_features <- names(sort(cv, decreasing=TRUE))[1:top_n]
      mat <- mat[,base::eval(top_cv_features)]
      return(mat)
    } else {
      return(mat) # no extra computation needed
    }
}
#' @name getTopFeatures
#' @rdname getTopFeatures
#' @title Select top features
#' @description Given a results data.frame, select the top features and return the up regulated, down regulated, and all features.
#' @param \code{results} Result of the formatLMMResults function. A data.frame that is formatted with required column names: "Feature", "Estimate", "P", \code{p_adjust_column}
#' @param \code{n_features} Number of features to be selected for the up-regulated or down-regulated selection. Default: 10
#' @param \code{est_thr} Threshold for Estimate. Default: 0
#' @param \code{p_adjust_column} Adjusted P-value column name in the results data.frame or the selected p.adjust method. Default: "FDR"
#' @param \code{p_adjust_thr} Threshold for the adjusted P-value method. Default: 0.05
#' @details The function intakes the results dataframe of the \code{formatLMMResults} function, performed after a \code{mixedModelIDE}, to return a list containing the top up-regulated, down-regulated, and all top features.
#' @return A list containing three character vectors: up, down, all
#' @export
#' @examples
#' library(GeomxTools)
#' library(sDAS)
#' data("kidney_targets")
#' geneSetList <- list(Keratins = c('KRT4','KRT5','KRT8','KRT9','KRT18','KRT19'),
#'                    CellMarkers = c('CD3E','CD68','CD4','PECAM1','PTPRC','FAP'),
#'                    "T-box" = c('TBX19', 'TBX2', 'TBX6', 'TBX3', 'EOMES'),
#'                    STAT = c('STAT3', 'STAT6', 'STAT2', 'STAT1', 'STAT5A', 'STAT5B'),
#'                    "Vitamin B6 metabolism" = c('PNPO','PSAT1','PHOSPHO2','PDXP','AOX1','PDXK'),
#'                    "Tyrosine metabolism" = c('COMT','FAH','ADH1A','MAOB','ALDH3B2','HGD','ALDH1A3','IL4I1','GOT1','ADH5','GOT2','MIF','ALDH3B1','AOX1','FAHD1','MAOA','HPD','DDC'),
#'                    "Nicotinate and nicotinamide metabolism" = c('NUDT12','NT5C','NAMPT','PNP','QPRT','NADK2','NMNAT1','NAPRT','BST1','NADSYN1','NMRK2','NT5C2','NT5C3B','NNMT','NADK','NNT','NT5C3A','NT5E','AOX1'))
#' geneSet_data <- geneSetAnalysis(object = kidney_targets,
#'                                 elt = "log_q",
#'                                 geneSet = geneSetList,
#'                                 species = "Hs")
#' pData(geneSet_data)$slide <- as.factor(pData(geneSet_data)$`slide name`)
#' pData(geneSet_data)$region <- as.factor(pData(geneSet_data)$region)
#' geneSetDE <- mixedModelDE(geneSet_data[, rownames(pData(geneSet_data))],
#'                            modelFormula = ~ region + (1 + region|slide),
#'                            groupVar = "region",
#'                            nCores=parallel::detectCores()-1,
#'                            multiCore = FALSE)
#' geneSet_results <- formatLMMResults(geneSetDE)
#' geneSet_results_getTopFeatures <- getTopFeatures(geneSet_results)
#' geneSet_results_getTopFeatures


getTopFeatures <- function(results,
                           n_features = 10,
                           est_thr = 0,
                           p_adjust_column = "FDR",
                           p_adjust_thr = 0.05) {

 # ########
 # Utility
 # ########
 if(!(p_adjust_column == toupper(p_adjust_column))){
   p_adjust_column <- toupper(p_adjust_column)
   message("Your p_adjust_column value has been formatted to be uppercase.")
 }

 # ################
 # Check User Input
 # ################

 # Ensure the results input is a dataframe
 if(!inherits(results, "data.frame")){
  stop("Please check the input; input must be a data.frame.")
 }

 # Ensure that the required columns are present in the results dataframe
 if(!("Feature" %in% colnames(results) & "Estimate" %in% colnames(results) & "P" %in% colnames(results))){
  stop('Please check the input dataframe; input dataframe does not contain the named columns expected.  \n')
 }
 if(!any(colnames(results) == p_adjust_column)){
   stop(paste0('Please check the input dataframe; the input dataframe does not contain the selected p_adjust_column, ', p_adjust_column, '.'))
 }
 # Ensure the n_features is a positive integer greater than 1
 if(is.null(n_features)){
     n_features = 10
     warning("Your number of features was NULL; n_features will default to 10.")
 }
 if(!is.null(n_features)){
  if (!inherits(n_features, "numeric")){
   stop("Please check the number of features; n_features must be numeric.")
  }
  if (n_features %% 1 != 0){
   stop("Please check the number of features; n_features must be an integer.")
  }
  if (n_features < 1){
   stop("Please check the number of features; n_features must be >= 1.")
  }
 }

 # Ensure that the estimate threshold is a number, 0 or greater
 if(is.null(est_thr)){
   est_thr = 0
   warning("Your estimate threshold value was NULL; est_thr will default to 0.")
 }
 if(!is.null(est_thr)){
  if (!inherits(est_thr, "numeric")){
   stop("Please check the estimate threshold; est_thr must be numeric.")
  }
  if (est_thr < 0){
   stop("Please check the estimate threshold; est_thr must be >= 0.")
  }
 }

# Ensure that the adjusted P-value column name is a character

# Ensure that the adjusted P-value threshold is a number, 0 or greater
 if(is.null(p_adjust_thr)){
     p_adjust_thr = 0.05
     warning("Your adjusted P-value threshold was NULL; p_adjust_thr will default to 0.05.")
 }
 if(!is.null(p_adjust_thr)){
  if (!inherits(p_adjust_thr, "numeric")){
   stop("Please check the adjusted P-value threshold; p_adjust_thr must be numeric.")
  }
  if (p_adjust_thr < 0){
   stop("Please check the adjusted P-value threshold; p_adjust_thr must be >= 0.")
  }
 }

 # ########
 # Process
 # ########

 results$invert_P <- -log10(results$P) * sign(results$Estimate)

 # Select top up-regulated features
 results <- results[order(results$invert_P, decreasing = TRUE), ]
 top_features <-
  list(up = subset(results, results[p_adjust_column] <= p_adjust_thr & Estimate > est_thr)[1:n_features, "Feature"])
 top_features$up <- top_features$up[!is.na(top_features$up)]

 # Select top down-regulated features
 results <- results[order(results$invert_P, decreasing = FALSE), ]
 top_features$down = subset(results, results[p_adjust_column] <= p_adjust_thr & Estimate < -1*est_thr)[1:n_features, "Feature"]
 top_features$down <- top_features$down[!is.na(top_features$down)]

 # Select top up and down regulated features
 top_features$all = c(top_features$up, rev(top_features$down))

 # #######
 # Return
 # #######

 return(top_features)
}


setGeneric(name="getUMAP", signature = c("object"),
           function(object, ...){
             standardGeneric("getUMAP")
           })
#' @name getUMAP
#' @rdname getUMAP
#' @title Wrapper function for umap's umap function
#' @description Run umap with additional options.
#' @param \code{object} a NanoStringGeoMxSet object.
#' @param \code{elt} \code{character} specifying the expression \code{matrix}.
#' @param \code{use_top_cv} Logical. Use the targets with the high coefficient of variation? Default is FALSE.
#' @param \code{top_n} Numeric integer >=1. Can be NULL. If \code{use_top_cv} is TRUE, use this many targets. If \code{top_n}
#' greater than the number of features, no subsetting will be used.
#' @param \code{log2_transform} Logical. To log2 transform data or not. Default is FALSE.
#' @param \code{umap_config} \code{umap.config} of arguments to pass to the umap function. See details.
#' @param \code{the_prefix} \code{character} giving the prefix to use for column names. Default is "UMAP".
#' @details This function is a wrapper for running the umap function with a \code{NanoStringGeoMxSet}
#' object with custom and built-in arguments. For the built-in arguments,
#' \code{getUMAP} passes the \code{umap.config} from umap::umap.defaults.
#' The user can optionally pass a \code{umap.config}
#'  object with different values. Note that the names of the element are required to be the same as the umap::umap.defaults version.
#'  Also note that no additional error or class checking is performed on those elements.
#'
#'  The umap values for samples will be added as a new assayDataElement with the number of columns
#'  equal to the number specified in umap_config$n_components (one row for each sample in object). The new element name will begin
#'  with \code{the_prefix}
#'
#' @return A modified \code{NanoStringGeoMxSet} object. The umaps values for samples will be added to pData(object)
#' with the number of new columns equal to the number specified in umap_config$n_components.
#'
#' @importFrom Biobase pData
#' @importFrom Biobase fData
#' @importFrom umap umap
#' @import GeomxTools
#' @export
#' @examples
#' library(GeomxTools)
#' kidney_targets <- updateGeoMxSet(kidney_targets)
#' getUMAP(kidney_targets, elt="log_q")

setMethod(f="getUMAP",
          signature = "NanoStringGeoMxSet",
          definition = function(object,
                                elt,
                                use_top_cv=FALSE,
                                top_n=100,
                                log2_transform=FALSE,
                                umap_config=umap::umap.defaults,
                                the_prefix="UMAP") {

  ## Check the input

  # elt must be a non-zero length character and specifies an expression matrix
  if(!inherits(elt, "character")){
    stop("Please provide a character to elt.")
  }
  if(nchar(elt) == 0){
      stop("elt can not be a zero-length variable name.")
  }
  if(!inherits(assayDataElement(object, elt = elt), "matrix")){
    stop("elt was not found in object. Please check spelling.")
  }

  # use_top_cv must be logical.
  if(!is.logical(use_top_cv)){
    stop("use_top_cv must be logical.")
  }

  # top_n must be NULL or a positive numeric integer. If NULL, use_top_cv must be FALSE.
  if(is.null(top_n) & use_top_cv){
    stop("Please provide the number of targets to use when restricting to the top CV targets.")
  }
  if(!is.null(top_n)){
    if(!inherits(top_n, "numeric")){
      stop("top_n must be numeric.")
    }
    if(top_n %% 1 != 0){
      stop(paste0("The top_n given, ", top_n, " is not an integer."))
    }
    if(top_n <= 0){
      stop("top_n must be >= 1.")
    }
  }

  # log2_transform must be logical. If it is TRUE, no elements
  # in the expression matrix may be zero.
  if(!is.logical(log2_transform)){
    stop("log2_transform must be logical.")
  }
  if(log2_transform & any(assayDataElement(object, elt=elt)<=0)){
    stop("Cannot do log2 transformation on values less than or equal to zero.")
  }

  # umap_config must be a list with the names below.
  if(!inherits(umap_config, "umap.config")){
    stop("umap_config must be a umap.config. See details.")
  }
  valid_names <-  names(umap::umap.defaults)

  if(!all(valid_names %in% names(umap_config))){
    missing_names <- valid_names[!valid_names %in% names(umap_config)]
    stop(paste0("Elements in umap_config must include elements in umap.defaults. ",
                paste(missing_names, collapse=", "), " were not found."))
  }

  # the_prefix must be a character string
  if(!inherits(the_prefix, "character")){
    stop("the_prefix must be a character.")
  }

  # New column names added (with the form <the_prefix>_coords_1, ..., <the_prefix>_coords_<umap_config$dim>)
  # must not already be present in the phenotype data.
  new_coord_names <- paste0(the_prefix, "_coords_",
                            as.character(seq.int(from=1, to=umap_config$n_components, by=1)))
  if(any(new_coord_names %in% colnames(pData(object)))){
    stop("One or more coordinate names already exists. Please specify a different the_prefix value.")
  }

  # Analysis
  if(log2_transform){
    mat <- t(log2(assayDataElement(object, elt = elt)))
  } else {
    mat <- t(assayDataElement(object, elt = elt))
  }
  if(use_top_cv){
    mat <- getTopCV(mat, top_n)
  }

  umap <- umap::umap(mat, config=umap_config)
  umap_layout <- umap$layout

  # Adjust column names.

  colnames(umap_layout) <- paste0(the_prefix, "_coords_", 1:ncol(umap_layout))
  # Add columns of pca$ind$coord to phenotype data.
  pData(object) <-  cbind(
      as(pData(object), "data.frame"),
      as.data.frame(umap_layout))

  # return object
  return(object)

})
setGeneric(name="makeQCHistogram", signature=c("object"),
           function(object, ...){
             standardGeneric("makeQCHistogram")
           })

#' @name makeQCHistogram
#' @rdname makeQCHistogram
#' @title Make QC Histograms
#' @description Plot histogram of QC metrics
#' @param \code{object} a NanoStringGeoMxSet object
#' @param \code{annotation_col} a single character vector to indicate which QC Metric to plot. Examples:'Trimmed \%', 'Stitched \%', 'Aligned \%', 'Saturated \%', 'area', 'nuclei', 'NegGeoMean'.
#' @param \code{bins} A single, numeric value to set the number of histogram bins.
#' @param \code{fill_by} Optional. A single character vector indicating which element of pData to facet and color the histograms (e.g. 'Segment').
#' @param \code{xintercept} Optional. A single, numeric value to set the x-intercept.
#' @param \code{scale_trans} Optional. Default (NULL) uses a linear scale. Other valid transformations from \code{ggplot2::scale_x_continuous(trans="")} (e.g. 'log10') can be passed.
#'
#' @details This function plots a histogram of QC Metrics that are stored in an NanoStringGeoMxSet object.
#'
#' The function takes a NanoStringGeoMxSet object with QC metric calculated and returns a histogram for each level of pData(object) specified.
#'
#' @return A \code{ggplot} object.
#'
#' @import GeomxTools
#' @export
#' @examples
#' library(GeomxTools)
#' makeQCHistogram(object = kidney, annotation_col = "Trimmed", bins = 50, fill_by = "segment", xintercept = 60)

setMethod(f = "makeQCHistogram",
          signature = "NanoStringGeoMxSet",
          definition = function(object,
                         annotation_col=NULL,
                         bins=NULL,
                         fill_by=NULL,
                         xintercept=NULL,
                         scale_trans=NULL) {

  ## Check the input

    # Object must be of class NanoStringGeoMxSet.
    if(!inherits(object, "NanoStringGeoMxSet")){
      stop(paste0(object, "must be of class NanoStringGeoMxSet"))
    }

    # annotation_col must a character present in colnames(sData(object)).
    if(is.null(annotation_col)){
      stop("Please specify an annotation_col in colnames(sData(object)) to plot")
    } else if(!annotation_col %in% colnames(sData(object))){
      stop(paste0(annotation_col, " must be present in colnames(sData(object))"))
    }

    # annotation_col must a character present in colnames(sData(object)).
    if(!inherits(sData(object)[,annotation_col],"numeric")){
      stop("values within annotation_col must be numeric")
    }

    # bins must must be a single integer if not NULL
    if(!is.null(bins) && !inherits(bins, "numeric")){
      stop("bins needs to be a single integer.")
    } else if(!is.null(bins) && length(bins)!=1){
      stop("bins needs to be a single integer.")
    } else if(is.null(bins)){
      warning("Bins not specified. Using bins=50.")
      bins<-50
    }

    # fill_by must a character present in colnames(pData(object)).
    if(!is.null(fill_by) && !inherits(fill_by, "character")){
      stop("fill_by needs to be a character present in pData(object).")
    } else if(!is.null(fill_by) && !fill_by %in% colnames(pData(object))){
        stop(paste0("fill_by is not NULL but ", fill_by, " is not present in pData(object)."))
    }

    # xintercept must be a single integer if not NULL.
    if(!is.null(xintercept) && !inherits(xintercept, "numeric")){
      stop("xintercept needs to be a single numeric value.")
    } else if(!is.null(xintercept) && length(xintercept)!=1){
      stop("xintercept needs to be a single numeric value.")
    }

    # scale_trans must a character present in  ggplot2::scale_x_continuous(trans=).
    if(!is.null(scale_trans) && !scale_trans %in% c("asn", "atanh", "boxcox", "date", "exp", "hms",
                                                    "identity", "log", "log10", "log1p", "log2", "logit",
                                                    "modulus", "probability", "probit", "pseudo_log", "reciprocal",
                                                    "reverse", "sqrt","time")){
      stop("scale_trans must be  option from ggplot2::scale_x_continuous(trans=)")
    }

  ## Plotting
  # Extract data frame of QC data from the NanoStringGeoMxSet object.
  plot_df <- sData(object)

  # Create histogram
  plt <- ggplot(plot_df,
                aes_string(x = paste0("unlist(`", annotation_col, "`)"),
                           fill = fill_by)) +
    geom_histogram(bins = bins) +
    geom_vline(xintercept = xintercept, lty = "dashed", color = "black") +
    theme_bw() + guides(fill = "none") +
    labs(x = annotation_col, y = "Segments, #", title = annotation_col)

  # Facet the histogram if "fill_by" is specified
  if(!is.null(fill_by)) {
    plt <- plt +
      facet_wrap(as.formula(paste("~", fill_by)), nrow = length(unique(plot_df[,fill_by])))
  }

  # Add continuous x-axis if "scale_trans" is specified
  if(!is.null(scale_trans)) {
    plt <- plt +
      scale_x_continuous(trans = scale_trans)
  }

  # Plot Histogram
  return(plt)

})
#' @name makeVolcano
#' @rdname makeVolcano
#' @title Generate Volcano Plot
#' @param df is a data.frame or tibble with columns
#'  Comparison, Feature, Estimate, and P. See details.
#' @param p_adjust_column character. Column name of the P-adjust data in \code{df}. For example, 'FDR'. Default is FDR.
#' @param to_label a character string of features to label.
#' @param label_color character The color for features of interest. Default = NULL.
#' @param SCALE numeric. The scaling factor. Used to extend the x-axis.
#' @param LWD numeric. The line width.
#' @param the_name character. A name to add to the plot. If NULL, will plot the Comparison name.
#' @param log_type a character for plotting the x-axis. Either 'log2(FC)' or 'FC'.
#' @param fc_cutoff numeric. The estimate cutoff (vertical lines)
#' @param pval_cutoff numeric. The horizontal line.
#' @details In addition to the four core columns needed in \code{df} (Contrast, Feature, Estimate,
#' P), an additional column is needed for the P-value adjustment if \code{p_adjust_column} is not "NONE".
#' Parameter \code{p_adjust_column}
#' references the column in \code{df} that is to be used (e.g., 'FDR', 'BY'). This column can not be 'P'.
#' @return A list containing the ggplot object and the data.frame that was used to generate the plot.
#' @importFrom ggrepel geom_text_repel
#' @importFrom dplyr filter
#' @export
#' @examples
#' library(plyr)
#' library(GeomxTools)
#' data("kidney_targets")
#' ## Generate DE results from a small set of genes
#' kidney_targets_subset <- kidney_targets[1:25, ]
#' kidney_targets_subset$slide <- factor(kidney_targets_subset$slide)
#' kidney_targets_subset$class <- factor(kidney_targets_subset$class)
#' de <- mixedModelDE(kidney_targets_subset,
#'                   elt = "log_q",
#'                   modelFormula = ~ class + (1|slide),
#'                   groupVar = "class",
#'                   nCores = parallel::detectCores()-1,
#'                   multiCore = FALSE)
#' de_format <- formatLMMResults(de)
#' # Volcano plot
#' makeVolcano(de_format, to_label=de_format$Feature[1:2],
#' p_adjust_column = 'FDR', label_color="grey",
#'   SCALE=1.4, LWD=1, the_name = NULL,
#'   log_type="log2(FC)", fc_cutoff=0.5, pval_cutoff=0.05)

makeVolcano <- function(df,
                        to_label=NULL,
                        p_adjust_column = "FDR",
                        label_color=NULL,
                        SCALE=1.4, LWD=1,
                        the_name = NULL,
                        log_type='log2(FC)',
                        fc_cutoff=1,
                        pval_cutoff=0.05){

 # ###################
 # User input checking
 # ###################

 # df needs to inherit data.frame
 if(!inherits(df, "data.frame")){
  stop("df needs to be a data.frame.")
 }
 # Check df columns
 if(!all(c('Feature', 'Comparison', 'Estimate', "P") %in% colnames(df))){
  stop("Check column names in df.")
 }
 # The p_adjust_column must be present in df and must not be 'P'.
 if(!any(colnames(df) == p_adjust_column)){
  stop(paste0("The column name provided, ", p_adjust_column, ", was not found in df."))
 }
 if(p_adjust_column=="P"){
  stop("p_adjust_column cannot be P since P is already used in the legend.")
 }
 # df needs observations
 if(nrow(df)==0){
  stop("No data found in df.")
 }
 # df cannot contain more than one comparison
 if(length(unique(df$Comparison))!=1){
  stop("Only one comparison can be plotted. Subset df.")
 }

 # label_color needs to be a single character value if to_label is not NULL.
 # Additionally, all values in to_label must be characters and present in df$Feature
 if(!is.null(to_label)){
  if(!inherits(label_color, "character")){
   stop("Please provide a character specifying the label color.")
  }
  if(length(label_color)!=1){
   stop("Please provide one label color for feature names.")
  }
  if(!inherits(to_label, "character")){
   stop("to_label must be a character vector.")
  }
  if(!all(to_label %in% df$Feature)){
   stop("Not all features in to_label are present in df. Check spelling.")
  }
 }

 # the_name must be a character of length one if not NULL
 if(!is.null(the_name)){
  if(!inherits(the_name, "character")){
   stop("the_name must be a character.")
  }
  if(length(the_name)!=1){
   stop("the_name must be a single value.")
  }
 }

 # log_type must be one of these two
 if(!log_type %in% c("FC", "log2(FC)")){
   stop("log_type needs to be FC or log2(FC)")
 }

 # SCALE, LWD, fc_cutoff, and pvalue_cutoff must all be a single, positive numeric
 if(!inherits(SCALE, "numeric")){
  stop("SCALE must be numeric.")
 }
 if(length(SCALE)!=1){
  stop("SCALE must be a single value.")
 }
 if(SCALE<0){
  stop("SCALE must be postive.")
 }
 if(!inherits(LWD, "numeric")){
  stop("LWD must be numeric.")
 }
 if(length(LWD)!=1){
  stop("LWD must be a single value.")
 }
 if(LWD<0){
  stop("LWD must be postive.")
 }
 if(!inherits(fc_cutoff, "numeric")){
  stop("fc_cutoff must be numeric.")
 }
 if(length(fc_cutoff)!=1){
  stop("fc_cutoff must be a single value.")
 }
 if(fc_cutoff<0){
  stop("fc_cutoff must be postive.")
 }
 if(!inherits(pval_cutoff, "numeric")){
  stop("pval_cutoff must be numeric.")
 }
 if(length(pval_cutoff)!=1){
  stop("pval_cutoff must be a single value.")
 }
 if(pval_cutoff<0){
  stop("pval_cutoff must be postive.")
 }

 # ##########
 # Processing
 # ##########

 dfx <- df
 # Parse the test name
 strings <- strsplit(dfx$Comparison[1], split=" vs ")[[1]]
 first <- strings[1]
 second <- strings[2]

 # Calculate the p_adjust_class
 p_adjust_class <- classifyPAdjustMethod(p_adjust_column)

 # Set significance levels
 dfx$sig_level <- paste0('NS or Estimate < ', fc_cutoff)
 dfx$sig_level[dfx$P < pval_cutoff & abs(dfx$Estimate)>=fc_cutoff] <- paste0('P < ', pval_cutoff)
 if(p_adjust_column == "NONE"){
  dfx$sig_level <- factor(dfx$sig_level, levels = c(paste0('NS or Estimate < ', fc_cutoff),
                                                    paste0('P < ', pval_cutoff)))
 } else {
  dfx$sig_level[dfx[,p_adjust_column] < 0.05 & abs(dfx$Estimate)>=fc_cutoff] <- paste0(p_adjust_class, ' < 0.05')
  dfx$sig_level[dfx[,p_adjust_column] < 0.001 & abs(dfx$Estimate)>=fc_cutoff] <- paste0(p_adjust_class, ' < 0.001')
  dfx$sig_level <- factor(dfx$sig_level, levels = c(paste0('NS or Estimate < ', fc_cutoff),
                                                    paste0('P < ', pval_cutoff),
                                                    paste0(p_adjust_class, ' < 0.05'),
                                                    paste0(p_adjust_class, ' < 0.001')))
 }

 # Set the significance colors
 the_colors <- c("orange2", "grey")
 names(the_colors) <- c(paste0('P < ', pval_cutoff), paste0('NS or Estimate < ', fc_cutoff))
 if(p_adjust_column != "NONE"){
  additional_colors <- c('#00708b', '#A6CE39')
  names(additional_colors) <- c(paste0(p_adjust_class, ' < 0.001'), paste0(p_adjust_class, ' < 0.05'))
  the_colors <- c(additional_colors, the_colors)
 }

 # Set the x-axis name
 if(log_type=="log2(FC)"){
   the_x_lab <- substitute(a %<-% s %->% b, list(a = paste0(first, "  "), b = paste0("  ", second), s = bquote("log"[2]*"FC")))
 } else {
    the_x_lab <- substitute(a %<-% " FC " %->% b,
         list(a = paste0(first, "  "), b = paste0("  ", second)))
 }

 # ########
 # Plotting
 # ########

 plt <- ggplot(data=dfx, aes(x=Estimate, y=-log10(P), label=Feature)) +
   geom_point(aes(colour=sig_level), alpha=0.75, size = 1.5) +
   scale_colour_manual(values=the_colors) +
   guides(colour=guide_legend(title="Significance\nLevel"))
 plt <- plt +
   geom_vline(xintercept=c(-fc_cutoff, fc_cutoff), lty = "dashed", lwd = LWD, color = "gray") +
   geom_hline(yintercept = -log10(pval_cutoff), lty = "dashed", lwd = LWD, color = "grey") +
   theme_bw(base_size = 14) +
   ylab(expression ('Significance, -log'[10]*'(P-value)')) +
   xlab(the_x_lab) +
   scale_x_continuous(limits=c(min(dfx$Estimate)*SCALE, max(dfx$Estimate)*SCALE)) +
   if(is.null(the_name)){
    ggtitle(label = paste0(dfx$Comparison[1]))
   } else{
    ggtitle(label = the_name)
   }

  if(!is.null(to_label)){
      plt <- plt + ggrepel::geom_text_repel(
      data=dfx %>% dplyr::filter(Feature %in% to_label),
        aes(x=Estimate, y=-log10(P), label=Feature),
      color = label_color, box.padding = 0.6, point.padding = 0.2,
      min.segment.length = 0.25, fontface = "bold", max.overlaps = 40)
  }

  # Add title
  if(!is.null(the_name)){
     plt <- plt + ggtitle(the_name)
  }

  return(list(plt, dfx))
}
setGeneric(name="performLinearModelDE", signature=c("object"),
           function(object, ...){
             standardGeneric("performLinearModelDE")
           })

#' @name performLinearModelDE
#' @rdname performLinearModelDE
#' @title Differential expression using a linear model
#' @description Perform differential expression on NanoString GeoMx feature level count data with a linear model.
#' @param \code{object} a NanoStringGeoMxSet object
#' @param \code{feature} a character vector indicating which data.frame in \code{assayData(object)} to perform differential expression on.
#'  rownames should indicate features (e.g. genes) and colnames should indicate sample ids
#' @param \code{annotations} a data.frame within \code{object} that includes sample annotations
#' @param \code{fixed_effect} the column name in \code{annotations} to be used as the fixed effect
#' @param \code{n_processors} the number of processors to use.
#' @param \code{combos} a list of fixed effect combinations to perform pairwise contrasts (e.g., list(c("Condition1", "Condition2"), c("WildType", "Mutant")))
#' @param \code{base_level} a vector of reference levels for each combo, listed in the same order as each combo. If NULL, the first value of the combo is used as the base level.
#' @param \code{formula} a character string for the linear model formula used.
#' @param \code{pAdjust} Optional. A character specifying which multiple testing correction method to use. Options are those from \code{p.adjust}
#' @param \code{log2_transform} logical. Whether to log2 transform the data (default) or not (FALSE).
#'
#'
#' @details This function performs differential expression using feature level counts stored in a NanoStringGeoMxSet object. This function is intended
#' for use cases where a random effect is not necessary, e.g. tissue microarrays (TMA). For more complex experimental designs, where the use of a
#' random effect is warranted, please use \code{\link[GeomxTools]{GeomxTools::mixedModelDE}} instead.
#'
#' @return A list of results for each set of fixed effect combinations.
#'
#' @import GeomxTools
#' @import parallel
#' @export
#' @examples
#' library(GeomxTools)
#' library(sDAS)
#' data("kidney_targets")
#' single_slide <- pData(kidney_targets)$`slide name` == "normal3"
#' single_slide <- kidney_targets[, single_slide]
#' lmResults <- performLinearModelDE(object = single_slide,
#'              feature = "exprs",
#'              annotations = pData(single_slide),
#'              fixed_effect = "region",
#'              combos = list(c("glomerulus", "tubule")),
#'              base_level = c("tubule"),
#'              n_processors=1,
#'              formula = ~ region,
#'              pAdjust = "BH",
#'              log2_transform = TRUE)
#' lmResults <- lmResults[[1]]

setMethod(f = "performLinearModelDE",
          signature = "NanoStringGeoMxSet",
          definition = function(object = object,
                                feature = "exprs",
                                annotations = pData(object),
                                fixed_effect = NULL,
                                combos = NULL,
                                base_level = NULL,
                                n_processors = 1,
                                formula = NULL,
                                pAdjust = "BH",
                                log2_transform = FALSE) {

     ## Check the input

     # Object must be of class NanoStringGeoMxSet.
     if(!class(object)[1] == "NanoStringGeoMxSet"){
      stop("object must be of class NanoStringGeoMxSet")
     }

     # Annotations must be of class data.frame
     if(!class(annotations) == "data.frame"){
      stop("annotations must be of class data.frame")
     }

     # feature must be an option from assayDataElementNames(object)
     if(!feature %in% assayDataElementNames(object)){
      stop(paste0(feature, " must be present in assayDataElementNames(object)"))
     }

     # Sample names must match between annotations and expression data
     if(!all(row.names(annotations) == colnames(assayDataElement(object, feature)))){
      stop("Sample names in annotations do not match expression data.")
     }

     if(!all(colnames(assayDataElement(object, feature)) == row.names(annotations))){
      stop("Sample names in expression matrix do not match those in annotations data")
     }

     # Fixed effect should be in the annotations and should contain at least 2 unique levels.
     if(is.null(fixed_effect)){
      stop("Please specify a fixed effect")
     } else if(!fixed_effect %in% colnames(annotations)){
      stop("the fixed effect must be a column name in the sample annotations")
     }

     if(length(unique(annotations[,fixed_effect])) < 2){
      stop("The fixed effect must have at least 2 levels to compare.")
     }

     # Combos needs to be specified and contain levels of fixed_effect
     if(is.null(combos)){
      stop("Please specify combos to test.")
     } else if(!class(combos) == "list"){
      stop("combos must be provided as a list")
     } else if(!all(unlist(combos) %in% unique(annotations[,fixed_effect]))){
      stop("the fixed effect levels must be included in combos")
     }

     # Base_level must be a vector and the length of base_level needs to be equal to the number of combos
     if(!is.null(base_level)){
      if(!is.vector(base_level)){
       stop("Base_level input must be a vector.")
      } else if(length(combos) != length(base_level)){
         stop("Each combo must be associated with a base_level for comparison.")
        }
     }

     # Each base_level must be a value of its respective combo
     if(!is.null(base_level)) {
      for (i in 1:length(combos)){
       if(!(base_level[i] %in% combos[[i]])){
        stop("Please check your base_level and combos. Base_level must be a value of its combo.")
       }
      }
     }

     # n_processors should be a single integer and class should be numeric
     if(!inherits(n_processors, "numeric")){
      stop(paste0("n_processors must be numeric"))
     } else if(n_processors %% 1 != 0){
      stop(paste0("n_processors must be an integer."))
     }

     if(n_processors > parallel::detectCores()){
      stop(paste0(n_processors, " processors specified but only ", parallel::detectCores(), " are available"))
     }

     if(n_processors < 1L){
      stop(paste0("n_processors must be between 1-", parallel::detectCores()))
     }

     # Check if there is a formula provided and if it contains a random effect
     if(is.null(formula)){
      stop("Please specify a model formula")
     } else if(grepl("\\|", as.character(formula)[2])){
      stop("The formula includes a random effect. Please use GeomxTools::mixedModelDE() instead.")
     }

     # The formula needs to contain the fixed effect
     model_variables <- all.vars(formula)
     if ("1" %in% model_variables) {
      model_variables <- model_variables[which(!(model_variables %in% "1"))]
     }
     if (!fixed_effect %in% model_variables){
      stop ("fixed_effect needs to be specified in the formula.")
     }
     # All model variables must be in the annotations and contain at least two levels.
     if(!all(model_variables %in% colnames(annotations))){
      stop("all model variables must be column names in the sample annotations")
     }
     for(i in model_variables){
      if(length(unique(annotations[,i])) < 2){
       stop("The model variables must contain at least 2 levels.")
      }
     }

     # Valid p.adjust method
     if(!pAdjust %in% c("bonferroni", "holm", "hochberg", "hommel", "BH", "fdr", "BY", "none")){
      stop("Multiple testing correction method must be an option from p.adjust")
     }

     # Is log2_transform logical?
     if(!inherits(log2_transform, "logical")){
      stop(paste0("log2_transform must be logical, e.g. TRUE or FALSE"))
     }


     ## Get feature data
     feature <-assayDataElement(object, feature)

     ## Subset annotations and convert model variables to factors
     annotations <- annotations[, model_variables, drop = FALSE]
     for (i in names(annotations)){
      if (inherits(i, "character")) {
       annotations[, i] <- as.factor(annotations[, i])
      }
     }

     ## Run the linear model
     resultsLM <- lapply(1:length(combos), function(i){

      # Identify pairs for contrasts
      a_combo <- combos[[i]]
      if(length(a_combo)>2){
       stop(paste0("combos should specify 2 levels to compare. Combination #", i, " contains ", length(a_combo), " levels."))
      }
      print(a_combo)

      # Identify the base_level for the comparison
      if (!is.null(base_level)) {
       a_base_level <- base_level[i]
      } else {
       a_base_level <- a_combo[1]
      }


      # Filter down to the samples that are within a_combo
      to_keep <- which(annotations[,which(colnames(annotations)==fixed_effect)] %in% a_combo)
      annotations <- annotations[to_keep, , drop= FALSE]
      annotations$Sample_ID <- row.names(annotations)
      feature <- feature[,annotations$Sample_ID]

      # Define the features.
      the_features <- rownames(feature)

      # set formula for below
      the_formula <- formula(paste("a_feature", as.character(formula)[2], sep = " ~ "))

      # Set up cluster
      cl <- parallel::makeCluster(n_processors)
      parallel::clusterExport(cl=cl,
                              varlist=c("feature", "annotations", "the_features", "formula", "a_combo", "log2_transform"),
                              envir=environment())

      # Process genes in parallel
      inner_res <- parallel::parLapply(cl, the_features, function(a_feature){

       print(a_feature)

       if(isTRUE(log2_transform)){
        featuresTransformed <- as.data.frame(log2(feature[a_feature,]))
       } else if (isFALSE(log2_transform)){
        featuresTransformed <- as.data.frame(feature[a_feature,])
       }

       colnames(featuresTransformed)[1] <- "a_feature"

       featuresTransformed$Sample_ID <- row.names(featuresTransformed)

       rownames(featuresTransformed) <- NULL

       dat <- base::merge(featuresTransformed, annotations, by="Sample_ID")
       if(var(dat$a_feature)==0) return(NULL)

       dat[[fixed_effect]] <- relevel(dat[[fixed_effect]], ref = a_base_level)

       model <- lm(as.formula(the_formula), data=dat)

       # Output
       cf <- data.frame(coefficients(summary(model)))[2,]
       cf$feature <- a_feature
       cf <- cf[,c(5,1:4)]
       colnames(cf) <- c("Feature", "Estimate", "SE", "tval", "P")
       cf$`-log10_pval` <- -log10(cf$P)

       comparison <- paste0(a_base_level, " vs ", a_combo[a_combo != a_base_level])
       cf$Comparison <- comparison
       row.names(cf) <- NULL
       res_list <- list(comparison=cf)
       return(res_list)
      })
      parallel::stopCluster(cl)

      # Tidy up the main summary results
      inner_res_summary <- do.call(rbind, lapply(inner_res, "[[", 1L))

      # P-adjustment based on subsets of data faceted by Comparison
      inner_res_summary$padj <- p.adjust(inner_res_summary$P, method = pAdjust)

      inner_res_summary <- inner_res_summary[, c("Feature", "Comparison", "Estimate","P", "padj")]
      colnames(inner_res_summary)[colnames(inner_res_summary)=="padj"] <- toupper(pAdjust)

      return(inner_res_summary)

     })
     return(resultsLM)

})
#' @name plotPairs
#' @rdname plotPairs
#' @title Pairwise plots
#' @param dat a data.frame.
#' @param color_by a column name in \code{dat} that is used for color.
#' @param color_scale Optional. A named vector specifying the manual color scheme. Default is NULL.
#' @return an object of class 'gg'.
#' @importFrom GGally ggpairs
#' @importFrom ggplot2 aes_string
#' @importFrom dplyr %>%
#' @export
#' @examples
#' library(GGally)
#' library(dplyr)
#' library(ggplot2)
#' library(GeomxTools)
#' data("kidney_targets")
#' example_dat <- pData(kidney_targets)[, c("GeneDetectionRate", "q_norm_qFactors", "class")]
#' example_color_scale <- c("DKD"="#3A6CA1", "normal"="orange")
#' plotPairs(example_dat, "class", example_color_scale)

plotPairs <- function(dat, color_by, color_scale=NULL){

 # dat must inherit a data.frame
 if(!inherits(dat, "data.frame")){
  stop("dat must inherit data.frame.")
 }

 if(!inherits(color_by, "character")){
  stop('color_by needs to be a character.')
 }

 if(length(color_by)!=1){
  stop("Only one column may be reference in color_by.")
 }

 # color_by must be a column in dat
 if(!color_by %in% colnames(dat)){
  stop(paste0("The specified column name, ", color_by, ", is not a column in dat."))
 }

 # color_by must be discrete
 if(!inherits(dat[,color_by], "character") & !inherits(dat[,color_by], "factor") & !inherits(dat[,color_by], "logical")){
  stop("color_by must be a character, factor, or logical data type.")
 }

 # Each unique value of color_by needs to be present in color_scale
 # if color_scale is not NULL
 if(!is.null(color_scale)){
  if(!inherits(color_scale, "character")){
   stop("color_scale must be NULL or a named character vector.")
  }
  if(!all(unique(dat[,color_by]) %in% names(color_scale))){
   stop(paste0("All unique values in ", color_by, " must be given a color."))
  }
 }

 n <- ncol(dat)
 p <- dat %>% ggpairs(.,
    mapping = aes_string(colour=color_by, alpha=0.5),
    columns=1:n, progress=FALSE,
    lower = list(continuous = wrap("smooth", alpha=0.3, size=0.3),
                 combo=wrap("facethist", bins=30))
  ) + theme_bw()

 if(!is.null(color_scale)){
  p <- p +
   scale_color_manual(values=color_scale) +
   scale_fill_manual(values=color_scale)
 }

 return(p)
}
#' @name plotStackedBars
#' @rdname plotStackedBars
#' @title plotStackedBars
#' @description plots stacked barplots of cell deconvolution data
#' @param object a NanoStringGeoMxSet object with feature type 'cellTypes'
#' @param elt the matrix to use for plotting. E.g., "beta"
#' @param label_name character. The name to use for the label. Default is NULL, which uses the name of \code{elt} value.
#' @param de_results optional data.frame from the function \code{\link{formatLMMResults}}.
#' @param to_label Either numeric or a character string of features to label. If numeric & de_results is not null, filters down to the to_label greatest fold change features.
#' @param ann_columns character vector giving the annotations in pData to use for labeling. Must have at least one annotation. Must not contain a column named "Sample_ID_Stacked_Bars".
#' @param other_name character giving the "other" category name. Default is "Other".
#' @param fill_colors a named character vector given the color for each cell type. Default (NULL) uses pals::alphabet
#' @param ann_colors a named list with colors. e.g., list('region'=c('one'='blue', 'two'='yellow'))
#' @param nuclei_column character. The column name in pData(object) that specifies nuclei. If present, plots the nuclei counts too. Default is NULL.
#' @seealso \code{\link{convertCellDecon}}
#' @importFrom pals alphabet
#' @examples
#' library(GeomxTools)
#' library(sDAS)
#' library(plyr)
#' library(dplyr)
#' data("kidney_cell_decon")
#' data("kidney_profile_matrix")
#' cd <- convertCellDecon(kidney_cell_decon, kidney_profile_matrix)
#' # format
#' cd <- cd[apply(assayDataElement(cd, elt="beta"), 1, sum)!=0,]
#' pData(cd)$region <- factor(pData(cd)$region)
#' pData(cd)$slide <- factor(pData(cd)$slide)
#' an_colors <- getColorPalette(pData(cd)[,c('region', 'slide')])
#' # color 4 cell types and have the rest be grey
#' to_label <- c("Glomerular.endothelium", "Podocyte",
#'               "Epithelial.progenitor.cell", "Connecting.tubule")
#' pal_cell_decon <- c("Glomerular.endothelium"="#A66293", "Podocyte"="#A66293",
#'               "Epithelial.progenitor.cell"="#D86769", "Connecting.tubule"="green", "Other"="grey")
#' plotStackedBars(cd, ann_columns = c("region", "slide"),
#'           to_label = to_label,
#'           fill_colors = pal_cell_decon)
#'
#' plotStackedBars(cd, ann_columns = c("region", "slide"),
#'           to_label = to_label,
#'           ann_colors=an_colors,
#'           fill_colors = pal_cell_decon)
#'
#' plotStackedBars(cd, ann_columns = c("region"),
#'           to_label = to_label,
#'           ann_colors=an_colors,
#'           fill_colors = pal_cell_decon)
#'
#' plotStackedBars(cd, ann_columns = c("region"),
#'           to_label = to_label,
#'           ann_colors=an_colors,
#'           fill_colors = pal_cell_decon,
#'           nuclei_column = "nuclei")
#'
#' plotStackedBars(cd, ann_columns = c("region"),
#'           elt = "prop_of_all",
#'           label_name = "Proportion of Cell Types",
#'           to_label = to_label,
#'           ann_colors=an_colors,
#'           fill_colors = pal_cell_decon,
#'           nuclei_column = "nuclei")
#' # Not Run
#' # # With differential abundance results
#' # da <- mixedModelDE(
#' #       object=cd,
#' #       elt="beta",
#' #       modelFormula = ~ region + (1+region|slide),
#' #       groupVar="region",
#' #       nCores=parallel::detectCores()-1,
#' #       multiCore=FALSE
#' #     )
#' # da <- formatLMMResults(da, p_adjust_method="BY")
#' #
#' # # The 5 cell types with the greatest fold change can be picked when
#' # # setting to_label to numeric and providing da results
#' # arrange(da, -abs(Estimate))$Feature[1:5]
#' # plotStackedBars(cd, ann_columns = c("region", "slide"),
#' #             elt = "prop_of_all",
#' #             label_name = "Proportion of Cell Types",
#' #             de_results = da,
#' #             to_label = 5,
#' #             ann_colors=an_colors,
#' #             nuclei_column = "nuclei")
#' @export

plotStackedBars <- function(object,
                      elt = "beta",
                      label_name = NULL,
                      de_results=NULL,
                      to_label=NULL,
                      ann_columns,
                      other_name="Other",
                      fill_colors=NULL,
                      ann_colors=NULL,
                      nuclei_column=NULL){

  if(!inherits(object, "NanoStringGeoMxSet")){
    stop('object needs be a NanoStringGeoMxSet.')
  }

  if(is.null(assayDataElement(object, elt=elt))){
    stop(paste0(elt, " is not an element in object."))
  }

  if(featureType(object)!="CellTypes"){
    stop(paste0("The feature type of object, ", featureType(object),
                ", is not CellTypes. Did you run convertCellDecon?"))
  }

  if(!is.null(label_name)){
    if(!inherits(label_name, "character")){
      warning("label_name will be converted to character.")
      label_name <- as.character(label_name)
    }
  } else {
    label_name <- elt
  }

  if(!is.null(de_results) & !is.null(to_label)){
    if(!inherits(to_label, "numeric")){
      stop("de_results is not NULL. Expecting a numeric for to_label")
    } else {
      foi <- arrange(as.data.frame(de_results), -abs(Estimate))[1:to_label, "Feature"]
    }
  } else if(inherits(to_label, "character")){
    foi <- to_label
  } else {
    foi <- rownames(object) # all features
  }


  if(is.null(ann_columns)){
    stop("ann_columns cannot be NULL. Please provide at least one annotation.")
  } else if(!all(ann_columns %in% colnames(pData(object)))){
      stop("ann_columns must refer to annotations in pData(object).")
  } else if("Sample_ID_Stacked_Bars" %in% ann_columns){
    stop("Sample_ID_Stacked_Bars cannot be a ann_column.")
  }


  # If other_name is already a cell type and a foi, stop since
  # this is a special category refered to downstream
  if(other_name %in% rownames(object)){
    stop(paste0(other_name, " is already a cell type. Please choose a different name for other_name."))
  }


  # Pivot data
  exprs <- t(assayDataElement(object, elt=elt))
  exprs <- as.data.frame(exprs)
  exprs$Sample_ID_Stacked_Bars <- row.names(exprs)
  row.names(exprs) <- NULL
  exprs <- exprs %>% tidyr::pivot_longer(cols=!Sample_ID_Stacked_Bars, names_to="cell_type", values_to="value")
  exprs$cell_type[!exprs$cell_type %in% foi] <- other_name
  exprs <- ddply(exprs, .(Sample_ID_Stacked_Bars, cell_type), summarize, SumValue=sum(value))

  # Get grand totals:
  grand_totals <- arrange(ddply(exprs, .(cell_type), summarize, Total=sum(SumValue)), Total)

  # Factor levels conditional on whether other_name is
  # present.
  if(other_name %in% exprs$cell_type){
    factor_levels <- c(other_name, setdiff(grand_totals$cell_type, other_name))
  } else {
    factor_levels <- grand_totals$cell_type
  }
  exprs$cell_type <- factor(exprs$cell_type, levels=factor_levels, order=TRUE)

  # Add annotation columns to data.frame and arrange
  annots <- pData(object) %>% dplyr::select(eval(ann_columns))
  annots$Sample_ID_Stacked_Bars <- row.names(annots)
  row.names(annots) <- NULL
  exprs <- base::merge(exprs, annots, by="Sample_ID_Stacked_Bars")

  # make each annotation a factor and sort to get the correct order for Sample_ID_Stacked_Bars.
  # If f is already an ordered factor, leave it.
  for(f in ann_columns){
    if(!is.ordered(exprs[,f])){
      exprs[,f] <- factor(exprs[,f], levels=unique(exprs[,f]), order=TRUE)
    }
  }
  exprs <- exprs[do.call("order", exprs[,ann_columns, drop=FALSE]), ]

  # Order Sample_ID_Stacked_Bars and sort data based on this ordering
  sample_order <- (exprs %>% dplyr::select(Sample_ID_Stacked_Bars, eval(ann_columns)) %>% distinct())$Sample_ID_Stacked_Bars
  exprs$Sample_ID_Stacked_Bars <- factor(exprs$Sample_ID_Stacked_Bars, levels=sample_order, order=TRUE)
  exprs <- exprs[do.call("order", exprs[,c(ann_columns, "Sample_ID_Stacked_Bars"), drop=FALSE]), ]

  # browser()
  # Plot results
  if(is.null(fill_colors)){
    fill_colors <- pals::alphabet(length(factor_levels))
    names(fill_colors) <- factor_levels
  } else if(length(setdiff(names(fill_colors), other_name))!=length(foi)){
    stop("The number of custom colors must match the number of cell types plotted.")
  }

  the_color="black"
  p_main <- ggplot() +
    geom_bar(data=exprs,
             aes(x=Sample_ID_Stacked_Bars, y=SumValue, fill=cell_type),
             position = "stack", stat="identity", width=1) +
    ylab(label_name) +
    scale_fill_manual(values=fill_colors) +
  theme_bw() +
   theme(axis.text.x = element_text(color = the_color, angle = 0, hjust = 0.5, vjust = 1, face = "plain"),
        axis.title.x = element_text(color = the_color, angle = 0, hjust = 0.5, vjust = 1, face = "plain")) +
   scale_y_continuous(expand = expansion(mult = 0)) +
  guides(fill=guide_legend(ncol=1)) +
  theme(legend.text=element_text(colour=the_color),
        legend.title=element_text(colour=the_color),
        legend.background = element_rect(fill='transparent'), #transparent legend bg
      legend.box.background = element_rect(fill='transparent')) + #transparent legend panel)
  coord_flip()

  p_main_strip <- p_main +
    guides(fill="none") +
    theme(axis.title.y=element_blank(),
          axis.text.y=element_blank(),
          axis.ticks.y=element_blank()) +
    theme(
      panel.grid.major.y = element_blank(),
      panel.grid.minor.y = element_blank(),
      axis.text.x.bottom = element_text(angle = 0, hjust = 0, vjust = 0, color = the_color),
      panel.border = element_rect(fill = 'transparent', color = the_color),
      axis.ticks = element_line(color = the_color),
      axis.text = element_text(color = the_color),
      axis.title = element_text(color = the_color),
      panel.background = element_rect(fill = "transparent"), # bg of the panel
      plot.background = element_rect(fill = "transparent", color = NA), # bg of the plot
    )

  # iterate over the annotations specified and make new ggplot objects
  # Return a list lof lists (basic plot, plot with borders stripped)
  annots_plots_list <- lapply(ann_columns, function(i){
    exprs_i <- exprs
    exprs_i$X <- exprs_i[,eval(i)]

    p_i_basic <- ggplot() +
      geom_bar(data=exprs_i %>% dplyr::select(Sample_ID_Stacked_Bars, X) %>% distinct(),
               aes(x=Sample_ID_Stacked_Bars, fill=X, y=1), stat="identity", width=1) +
      guides(fill=guide_legend(title=i)) +
      theme_bw() +
      ylab(i) +
      theme(legend.text=element_text(colour=the_color),
        legend.title=element_text(colour=the_color),
        legend.background = element_rect(fill='transparent'), #transparent legend bg
      legend.box.background = element_rect(fill='transparent')) + #transparent legend panel)
      coord_flip()

    if(!is.null(ann_colors)){
      if(i %in% names(ann_colors)){
        p_i_basic <- p_i_basic +
          scale_fill_manual(values=ann_colors[[i]])
      }
    }

    p_i_strip <- p_i_basic +
      guides(fill="none") +
      theme(axis.title.y=element_blank(),
        axis.text.y=element_blank(),
        axis.ticks.y=element_blank()) +
      theme(#axis.title.x=element_blank(),
        axis.text.x=element_blank(),
        axis.ticks.x=element_blank()) +
      theme(
        panel.grid.major.y = element_blank(),
        panel.grid.minor.y = element_blank(),
        panel.border = element_rect(fill = NA, color = NA),
        axis.ticks = element_line(color = the_color),
        axis.text = element_text(color = the_color),
        axis.title = element_text(color = the_color),
        panel.background = element_rect(fill = "transparent"), # bg of the panel
        plot.background = element_rect(fill = "transparent", color = NA), # bg of the plot
      )

    return(list(p_i_basic, p_i_strip))
  })

  # Add nuclei plot, if applicable
  nuc_needs_plotting <- FALSE
  if(!is.null(nuclei_column)){
    if(!nuclei_column %in% colnames(pData(object))){
      stop(paste0("The nuclei column, ", nuclei_column, " was not found."))
    } else if(!inherits(pData(object)[,nuclei_column], "numeric")){
      stop(paste0("The nuclei_column is not numeric."))
    }

    nuc_needs_plotting <- TRUE

    nuc_df <- data.frame(Sample_ID_Stacked_Bars=colnames(object), nuclei=pData(object)[,nuclei_column])
    nuc_df$Sample_ID_Stacked_Bars <- factor(nuc_df$Sample_ID_Stacked_Bars, levels=levels(exprs$Sample_ID_Stacked_Bars), order=TRUE)
    nuc_df <- nuc_df %>% arrange(Sample_ID_Stacked_Bars)

    p_nuc_basic <- ggplot() +
      geom_point(data=nuc_df,
               aes(x=Sample_ID_Stacked_Bars, y=nuclei), stat="identity") +
      theme_bw() +
      ylab("Nuclei") +
      theme(legend.text=element_text(colour=the_color),
        legend.title=element_text(colour=the_color),
        legend.background = element_rect(fill='transparent'), #transparent legend bg
      legend.box.background = element_rect(fill='transparent')) + #transparent legend panel)
      coord_flip()

    p_nuc_strip <- p_nuc_basic +
      guides(fill="none") +
      theme(axis.title.y=element_blank(),
        axis.text.y=element_blank(),
        axis.ticks.y=element_blank()) +
      theme(
        panel.grid.major.y = element_blank(),
        panel.grid.minor.y = element_blank(),
        panel.border = element_rect(fill = 'transparent', color = the_color),
        axis.ticks = element_line(color = the_color),
        axis.text = element_text(color = the_color),
        axis.title = element_text(color = the_color),
        panel.background = element_rect(fill = "transparent"), # bg of the panel
        plot.background = element_rect(fill = "transparent", color = NA), # bg of the plot
      )

  }


  # Combine individual plots
  plot_list <- lapply(annots_plots_list, "[[", 2L)
  plot_list[['main']] <- p_main_strip
  if(nuc_needs_plotting){
    plot_list[['nuclei']] <- p_nuc_strip
  }

  plot_legends <- lapply(annots_plots_list, function(i){
    return(cowplot::get_legend(i[[1]] +
     theme(legend.text=element_text(colour="black"),
         legend.title=element_text(colour="black"))))
  })
  plot_legends[['main']] <- cowplot::get_legend(p_main + theme(legend.text=element_text(colour="black"),
         legend.title=element_text(colour="black")))

  if(nuc_needs_plotting){
    the_rel_widths <- c(rep(1, length(plot_list)-2), 6, 2)
  } else {
    the_rel_widths <- c(rep(1, length(plot_list)-1), 10)
  }
  p_first_row <- cowplot::plot_grid(
    plotlist=plot_list,
    rel_widths = the_rel_widths,
    align="vh", nrow=1
  )
  p_second_row <- cowplot::plot_grid(
    plotlist=plot_legends,
    rel_widths = c(rep(1, length(plot_legends))),
    nrow=1
  )

  p_combined <- cowplot::plot_grid(p_first_row, p_second_row, nrow=2, rel_heights=c(1.5,1))
  return(p_combined)


}







#' @description Function to format list of characters. c('a', 'b', 'c') => a, b, and c
formatList <- function(vec){
  if(length(vec)==1){
    return(vec)
  } else if(length(vec)==2){
    return(paste0(vec[1], " and ", vec[2]))
  } else {
    vec_minus_one <- vec[1:(length(vec)-1)]
    vec_last_one <- rev(vec)[1]
    return(paste0(paste0(vec_minus_one, collapse = ", "), ", and ", vec_last_one))
  }
}

#' @description Function for programmatically providing the features and 
#' pathways of interest.
callFOI <- function(foi, poi, featureType="genes"){
  if(!featureType %in% c("genes", "proteins")){
    stop("featureType must be genes or proteins.")
  }
  if(length(foi)>0 | length(poi)>0){
    msg <- "Additionally, we tracked "
    if(length(foi)>=10){
      msg <- paste0(msg, "the ", length(foi), " ", featureType)
    } else if(length(foi)>1){
      msg <- paste0(msg, formatList(foi), " ", featureType)
    } else if(length(foi)>0){
      msg <- paste0(msg, "the ", foi, " ", gsub('.{1}$', '', featureType))
    }
    if(length(foi)>0){
      msg <- paste0(msg, " of interest")
    }
    if(length(poi)==0){
      msg <- paste0(msg, " throughout the analysis.")
    } else {
      if(length(poi)>=2){
        msg <- paste0(msg, " and the ", length(poi), " pathways ")
      } else {
        msg <- paste0(msg, " and the ", poi, " pathway ")
      }
      msg <- paste0(msg, "of interest throughout the analysis.")
    }
    return(msg)
  } 
}

printFactors <- function(facts){
  if(length(facts)>0){
    msg <- paste0("Individual ROI/AOI segments are annotated by ", formatList(facts), ".")
    return(msg)
  }
}

getSuffix <- function(x){
  if(x<0 | x>1){
    stop("A number between 0 and 1 is need.")
  }
  x <- 100*x
  i <- x %% 10
  j <- x %% 100
  if(i == 1 && j!=11){
    return(paste0(x, "st"))
  } else if(i==2 && j!=12){
    return(paste0(x, "nd"))
  } else if(i==3 && j!=13){
    return(paste0(x, "rd"))
  } else {
    return(paste0(x, "th"))
  }
}

#' @description Conditional statement displaying the number of features or segments 
#' that are removed. Only prints; does not filter.
printRemoved <- function(names, type){
  if(!type %in% c("features", "segments")) stop("type must be features or segments.")
  if(length(names)>0){
    if(length(names)==1){
      return(paste0("One ", gsub(".{1}$", "", type), " was removed by the above QC."))
    } else {
      return(paste0("A total of ", length(names), " ", type, " were removed by the above QC."))
    }
  } else {
    return(paste0("No ", type, " were removed in the above section."))
  }
}

printCellProfileInfo <- function(mat, is_custom, url=NULL, author=NULL){
  if(!is_custom){
    if(is.null(url)){
      stop("No URL given.")
    }
    return(paste0("For this analysis, we used [", ifelse(is.null(author), "this profile matrix", author), "](", url, ") which contains ", ncol(mat), " cell types."))
  } else {
    return(paste0("For this analysis we integrated the custom profile matrix that was supplied that contains ", ncol(mat), " cell types."))
  }
}



get_snr <- function(exp, backgrounds, x=1, type="SBR", n_processors=4){
  # exp has features as columns
  # type needs to be:
  #      "SBR": signal (S) to background (B) ratio = signal / background: 
  #      "SSR": signal to standard deviation ratio (S-B)/sigma_b 
  # n_processors number of processors to use
  # backgrounds are the control molecules.
  # signals are the signal molecules.
  if(!type %in% c("SBR", "SSR")){
    stop("Needs SBR or SSR for type")
  }
  background_features <- exp[, backgrounds]
  features <- exp
  # We want to keep the backgrounds for comparison
  # to_keep <- setdiff(colnames(features), backgrounds)
  # features <- features[,to_keep]
  
  if(class(backgrounds)=="numeric"){
    background_geoMean <- background_features
    background_sd <- rep(0, length(background_features))
  } else {
    background_geoMean <- apply(background_features, 1, EnvStats::geoMean)
    background_sd <- apply(background_features, 1, sd)
    # Rare cases, if sd = 0, SSR will not work. Replace with the next lowest SD
    background_sd[background_sd==0] <- as.numeric(sort(background_sd[background_sd>0])[1])
  }
  cl <- parallel::makeCluster(n_processors)
  parallel::clusterExport(cl=cl, varlist=c("x", "type", "features", "type", "background_geoMean", 
                                           "background_sd"), envir=environment())
  
  inner_res <- parLapply(cl, 1:ncol(features), function(j){
    if(type=="SBR"){
      out <- as.numeric((features[,j] / background_geoMean))
    } else if(type=="SSR"){
      out <- as.numeric((signif(features[,j], 5) - signif(background_geoMean, 5))/(x*background_sd)) # sig. digits.
    } else {
      stop("something went wrong. error code 1.")
    }
    #snr[,j] <<- as.numeric((exp[,j] - sample_LOD) / sample_geo_means)
    out <- data.frame(out)
    colnames(out) <- colnames(features)[j]
    rownames(out) <- rownames(features)
    return(out)
  })
  
  stopCluster(cl)
  the_snr <- do.call(cbind, inner_res)
  
  if(all(colnames(the_snr) == colnames(features)) == FALSE){
    warning("Column names were not in the original order and will be adjusted.")
    the_snr <- the_snr %>% dplyr::select(eval(colnames(features)))
  }
  
  if(any(is.na(the_snr))){
    stop("Some proteins not computed. Check function.")
  } else {
    return(data.matrix(the_snr))
  }
}

pairs_custom <- function(dat, column_of_interest){
  require(rlang)
  n <- ncol(dat)
  p <- dat %>% ggpairs(., 
                       mapping = ggplot2::aes(colour=as.character({{column_of_interest}}), alpha=0.1),
                       columns=1:n, progress=FALSE,
                       lower = list(continuous = wrap("smooth", alpha=0.3, size=0.3),
                                    combo=wrap("facethist", bins=30))
  ) + theme_bw() 
  return(p)
}