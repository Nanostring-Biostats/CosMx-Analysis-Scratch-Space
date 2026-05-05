#' Rank based inverse normal transformation
#'
#' This function is copied and borrowed from `RNOmni` package.
#' Applies the rank-based inverse normal transform (INT) to a numeric vector.
#' The INT can be broken down into a two-step procedure.
#' In the first, the observations are transformed onto the probability scale using the
#' empirical cumulative distribution function (ECDF).
#' In the second, the observations are transformed onto the real line, as Z-scores, using the probit function.
#'
#' @param u Numeric vector.
#' @param k Offset. Defaults to (3/8), correspond to the Blom transform.
#'
#' @return Numeric vector of rank normalized measurements.
#'
#' @examples
#'
#' # Draw from chi-1 distribution
#' y <- rchisq(n = 1e3, df = 1)
#' # Rank normalize
#' z <- RankNorm(y)
#' # Plot density of transformed measurement
#' plot(density(z))
#'
#' @export
#'

RankNorm <- function(u, k = 0.375) {
  if (!is.vector(u)) {
    stop("A numeric vector is expected for u.")
  }
  if ((k < 0) || (k > 0.5)) {
    stop("Select the offset within the interval (0, 0.5).")
  }
  if (sum(is.na(u)) > 0) {
    stop("Please exclude observations with missing measurements.")
  }
  n <- length(u)
  r <- rank(u)
  out <- qnorm((r - k) / (n - 2 * k + 1))
  return(out)
}


#' Run DE analysis models on SMI data.
#'
#' @param assay_matrix counts or normalized assay matrix.
#'                     Raw counts recommended for count distributions s.a. negative binomial (i.e., family = nbinom2).
#'                     Normalized is recommended for gaussian models.
#' @param metadata data.table or data.frame of meta.data associated with assay.
#'                 Must contain groupVar column and any covariates in the `formula` argument used for DE.
#'                 If cell id is not included as a column in metadata, should correspond to rownames in the metadata data.frame object.
#' @param formula right hand side of formula used in DE, possibly containing random effects.
#'                Must contain groupVar fixed effect.
#'                It's recommended that raw count models with a log-link (i.e., negative binomial, poisson)
#'                , use a corresponding offset for library size in formula
#'                (i.e., ~de_variable + offset(log(totalcounts))), where totalcounts is a column in the metadata.
#' @param groupVar optional categorical variable for which you may specify the ordering of factor levels, this column must exist in meta.data slot
#' @param groupVar_levels optional factor ordering for groupVar.
#' @param neighborhood_counts Deprecated argument (will be removed in the future); optional object created by smiDE::measure_neighbor_expression_by_celltype, which contains lists of
#'                            the expression counts in neighboring cells by cell type compared to a cell type of reference
#'                            (see `help(measure_neighbor_expression_by_celltype)`).
#'                            May be used to control for confounding or 'bleed-over' in DE model (see examples)
#' @param spatial_model optional list for fitting models with spatially correlated random effects.
#'                      This list must have "name" (indicating the type of spatial random effect model), as a minimum requirement.
#'                      See documentation and examples for details.
#'                      \link{spatial_model}
#' @param pre_de_obj optional list object created by smiDE::pre_de, which contains a table of cell_adjacencies `cell_adjacency_dt` used for calculating neighbor-cell expression,
#'                   and, if created with smiDE::pre_de(adjacencies_only = FALSE), an `nblist` object with pre-computed neighbor-cell expressions.
#'                   (see `help(pre_de)` or `help(measure_neighbor_expression_by_celltype)` for details).
#'                   May be used to control for confounding or 'bleed-over' in DE model (see examples)
#' @param neighbor_expr_overlap_weight_colname optional argument used in the case of provided `pre_de_obj` created with `cell_adjacencies_only = TRUE`.
#'                                             Used for neighbor expression control variable, default is NULL (unweighted), but can also be specified by a weight column in the `cell_adjacency_dt`
#'                                             (a column 'weight' is pre-computed in the `cell_adjacency_dt` as `1/distance` between neighboring cells).
#'                                             Cells with higher weights (closer to the modeled cells) can be given more
#'                                             weight in computing the neighbor expression of a given gene.  Coincides with assumption that closer cells are more likely to overlap/cause segmentation errors with the modeled cells.
#' @param neighbor_expr_overlap_agg whether to use the mean "sum" (sum total) or average "mean" of neighbor cell expression as the control variable. Default is "sum".
#' @param neighbor_expr_cell_type_metadata_colname column to use for determining neighbor cell expression by cell type.  Represented in pre_de_obj$cell_adjacency_dt
#' @param neighbor_expr_totalcount_normalize Defaults to TRUE, neighbor expression covariate is normalized by the library size or 'totalcounts' across all genes for a given cell.
#'                                           y_norm{cell, gene} = y_counts{cell, gene} * mean(totalcounts{allcells}) / totalcounts{cell}
#' @param neighbor_expr_totalcount_scalefactor If NULL, and neighbor_expr_totalcount_normalize is set to TRUE, this will calculate the cell-specific scaling factors for normalization directly from `assay_matrix`.
#' @param nCores = 1, number of cores to use, set to 1 if running in serial mode
#' @param verbose logical; default TRUE. If TRUE, prints per-target model fitting progress messages.
#' @param family a character string naming a family function,
#'               See `family` for a generic discussion of families.
#'               Default is nbinom2 (with log link), specifying Negative binomial distribution: quadratic parameterization (Hardin & Hilbe 2007). V = mu*(1+mu/phi) = mu+mu^2/phi.
#' @param targets optional vector of targets to run DE for.  If not supplied, run models for all targets in assay.
#' @param cellid_colname column name in metadata corresponding to cell id.
#' @param ... further arguments to be passed to model fitting function.
#'            See example below.
#'
#' @return data.table of DE results, with one row per target, marginal means and their SE's of the DE groups
#'         , estimated fold change and p-values for non-zero difference (identity link) or ratio \eqn{\neq 1} (log link) between groups.
#'
#' @examples
#'
#' library(Seurat)
#' library(data.table); setDTthreads(1)
#' datadir<-system.file("extdata", package = "smiDE")
#' sem <- readRDS(paste0(datadir, "/small_nsclc.rds"))
#' metainfo <- data.table(sem@meta.data)
#' totalcount_scalefactors <- mean(metainfo[["totalcounts"]]) / metainfo[["totalcounts"]]
#' names(totalcount_scalefactors) <- colnames(sem)
#' sem <- Seurat::SetAssayData(
#'   sem,
#'   "data",
#'   sem[["RNA"]]@counts %*% Matrix::Diagonal(x = totalcount_scalefactors, names = colnames(sem))
#' )
#'
#' pre_de_obj <- pre_de(
#'   metadata = metainfo,
#'   cell_type_metadata_colname = "cell_type",
#'   split_neighbors_by_colname = "tissue",
#'   mm_radius = 0.05,
#'   sdimx_colname = "sdimx",
#'   sdimy_colname = "sdimy",
#'   verbose = TRUE
#' )
#'
#' fibroblast_and_macrophage_cells <- metainfo[cell_type %in% c("fibroblast", "macrophage"), cell_ID]
#' de_results <- smi_de(assay_matrix = sem[["RNA"]]@counts,
#'   metadata = metainfo[cell_ID %in% fibroblast_and_macrophage_cells],
#'   formula = ~RankNorm(otherct_expr) + niche + tissue + offset(log(totalcounts)),
#'   family = "nbinom2",
#'   targets = rownames(sem)[1:5],
#'   pre_de_obj = pre_de_obj,
#'   neighbor_expr_cell_type_metadata_colname = "cell_type",
#'   neighbor_expr_overlap_weight_colname = NULL,
#'   neighbor_expr_overlap_agg ="sum",
#'   neighbor_expr_totalcount_normalize = TRUE,
#'   neighbor_expr_totalcount_scalefactor = totalcount_scalefactors,
#' )
#'
#' results(de_results, "pairwise", variable = "niche", targets = rownames(sem)[1:2])
#' results(de_results, "one.vs.rest", variable = "tissue", targets = rownames(sem)[1:2])
#'
#' ## fit a mixed model with all cells, cell type as a fixed effect covariate,
#' ## tissue/sample id as a random effect covariate
#'
#' de_results <- smi_de(assay_matrix = sem[["RNA"]]@counts,
#'   metadata = metainfo,
#'   formula = ~RankNorm(otherct_expr) + cell_type + (1 | tissue) + offset(log(totalcounts)),
#'   family = "nbinom2",
#'   targets = rownames(sem)[1:5],
#'   pre_de_obj = pre_de_obj,
#'   neighbor_expr_cell_type_metadata_colname = "cell_type",
#'   neighbor_expr_overlap_weight_colname = NULL,
#'   neighbor_expr_overlap_agg ="sum",
#'   neighbor_expr_totalcount_normalize = TRUE,
#'   neighbor_expr_totalcount_scalefactor = totalcount_scalefactors,
#' )
#'
#' results(de_results, "pairwise", variable = "cell_type", targets = rownames(sem)[1:2])
#' results(de_results, "one.vs.rest", variable = "cell_type", targets = rownames(sem)[1:2])
#'
#' @export
#' @import data.table
#' @importFrom stats as.formula coef family formula logLik qnorm terms.formula
#' @importFrom methods as
#' @importFrom utils stack
smi_de <- function(
  assay_matrix,
  metadata,
  formula,
  neighborhood_counts = NULL,
  pre_de_obj = NULL,
  groupVar = NULL,
  groupVar_levels = NULL,
  nCores = 1,
  verbose = TRUE,
  family = "nbinom2",
  targets = NULL,
  cellid_colname = "cell_ID",
  spatial_model = NULL,
  neighbor_expr_overlap_weight_colname = NULL,
  neighbor_expr_overlap_agg = c("sum", "mean"),
  neighbor_expr_cell_type_metadata_colname = "cell_type",
  neighbor_expr_totalcount_normalize = TRUE,
  neighbor_expr_totalcount_scalefactor = NULL,
  ...
) {

  metainfo <- data.table::copy(as.data.table(metadata))
  spatial_args <- NULL
  if (!is.null(spatial_model)) {
    stopifnot("spatial_model argument should be NULL,\n or a list including the 'name' of spatial model and corresponding model fitting arguments (see vignette for details)"
               = "name" %in% names(spatial_model))

    xycols <- c(spatial_model[["x_coord_col"]], spatial_model[["y_coord_col"]])
    if (spatial_model[["name"]] == "GP_INLA") {
      if (!requireNamespace("INLA", quietly = TRUE)) {
        stop("'GP_INLA' is specified but INLA package is not installed. To run a spatial model using INLA 'GP_INLA', please install the INLA R package: https://www.r-inla.org/download-install")
      }
      if (spatial_model[["name"]] == "GP_INLA") {
        if (!all(c("A", "priors") %in% names(spatial_model))) {
          if (!"mesh" %in% names(spatial_model)) {
            if (is.null(xycols)) {
             stop(paste0("If using a GP_INLA spatial random effect model, and not specifying 'priors', and 'A' matrix\n, need to specify 'x_coord_col' and 'y_coord_col' to proceed.  
                          Please specify using spatial_args = list(name = 'GP_INLA', x_coord_col = 'xname', y_coord_col = 'yname'"
                          ))
            }
            if (!all(c(xycols) %in% colnames(metadata))) {
              msg <- paste0("If using a GP_INLA spatial random effect model, and not specifying 'priors', and 'A' matrix\n, need to specify 'x_coord_col' and 'y_coord_col' to proceed.  
                            Expected X/Y Columns: "
                                   , xycols[1]
                                   , " and "
                                   , xycols[2]
                                   , " not found. Please specify using spatial_args = list(name = 'GP_INLA', x_coord_col = 'xname', y_coord_col = 'yname'"
                          )
              stop(msg)
            }
          }
        }
      }
    }

    default_args <- switch(spatial_model[["name"]]
                           , "GP_Matern" = list(fixed = list(nu = 0.5)
                                               , x_coord_col = "sdimx"
                                               , y_coord_col = "sdimy"
                                               , k_prop_n = 0.2
                                               , k = NULL
                                               , split_neighbors_by_colname = "Run_Tissue_name"
                                               , spatial_random_effect = ~Matern(1 | sdimx_cluster + sdimy_cluster %in% Run_Tissue_name)
                                               )
                           , "GP_INLA" = {
                                         ### some defaults here computed from the data
                                         rangey <- diff(range(metainfo[[spatial_model[["y_coord_col"]]]]))
                                         rangex <- diff(range(metainfo[[spatial_model[["x_coord_col"]]]]))
                                         ncells <- nrow(metainfo)
                                         list(x_coord_col = spatial_model[["x_coord_col"]]
                                             , y_coord_col = spatial_model[["y_coord_col"]]
                                             , max.edge.mesh = mean(c(rangex, rangey)) / sqrt(ncells) * c(5, 10)
                                             , offset.mesh = c(-0.05, -0.1)
                                             , alpha.inlaprior = 2
                                             , prior.range.inlaprior = c(mean(c(rangex, rangey)), 0.5)
                                             , prior.sigma.inlaprior = c(0.5, 0.5)
                                             , constr.inlaprior = TRUE
                                             , spatial_random_effect = ~f(s, model = spde)
                                             , quantiles = c(0.025, 0.5, 0.975)
                                             )
                             }
                           )

    ### start with formal-default arguments for package model fitting function
    formalss <-  switch(spatial_model[["name"]]
                        , "GP_Matern" = names(formals(spaMM::fitme))
                        , "GP_INLA" = names(formals(INLA::inla))
                        )

    ### override any formal-defaults with defaults specified above
    override_args <- intersect(names(spatial_model)
                               , formalss)
    argnames <- formalss
    spatial_args <- vector(mode = "list", length = length(argnames))
    names(spatial_args) <- argnames
    for (arg in names(default_args)) {
      spatial_args[[arg]] <- default_args[[arg]]
    }

    ### override defaults with user-passed arguments
    for (arg in names(spatial_model)) {
      spatial_args[[arg]] <- spatial_model[[arg]]
    }
    spatial_args[c("formula", "data", "family")] <- NULL
    for (arg in names(spatial_args)) {
      if (is.null(spatial_args[[arg]])) spatial_args[[arg]] <- NULL
    }

    xycols <- c(spatial_args[["x_coord_col"]], spatial_args[["y_coord_col"]])
    if (spatial_model[["name"]] == "GP_INLA") {
      if (!all(c("A", "priors") %in% names(spatial_args))) {
        if (!"mesh" %in% names(spatial_args))
          if (!all(c(xycols) %in% colnames(metadata))) {
            msg <- paste0("If using a GP_INLA spatial random effect model, and not specifying 'priors', and 'A' matrix\n, need to specify 'x_coord_col' and 'y_coord_col' to proceed.  
                          Expected X/Y Columns: "
                                 , xycols[1]
                                 , " and "
                                 , xycols[2]
                                 , " not found. Please specify using spatial_args = list(name = 'GP_INLA', x_coord_col = 'xname', y_coord_col = 'yname'"
                        )
            stop(msg)
          }
          spatial_args[["mesh"]] <- INLA::inla.mesh.2d(loc = metainfo[, c(spatial_args[["x_coord_col"]], spatial_args[["y_coord_col"]]), with = FALSE]
                                                       , offset = spatial_args[["offset.mesh"]]
                                                       , max.edge = spatial_args[["max.edge.mesh"]]
                                                       )

          spatial_args[["priors"]] <- INLA::inla.spde2.pcmatern(spatial_args[["mesh"]]
                                              , alpha = spatial_args[["alpha.inlaprior"]]
                                              , prior.range = spatial_args[["prior.range.inlaprior"]]
                                              , prior.sigma = spatial_args[["prior.sigma.inlaprior"]]
                                              , constr = spatial_args[["constr.inlaprior"]]
          ) # Need to specify priors (not sure what are the best default values)

          spatial_args[["A"]] <- INLA::inla.spde.make.A(mesh = spatial_args[["mesh"]]
                                                         , loc = data.matrix(metainfo[, xycols, with = FALSE]))
      }
      stopifnot(all(c("A", "priors") %in% names(spatial_args)))
    }

    if (spatial_args[["name"]] == "GP_Matern") {
      spamm_terms <- spatial_args[["spatial_random_effect"]]
      if (!all(all.vars(spamm_terms) %in% names(metainfo))) {
        if (!all(c(xycols) %in% colnames(metadata))) {
          msg <- paste0("Expected X/Y Columns: "
                               , xycols[1]
                               , " and "
                               , xycols[2]
                               , " not found. Please specify using spatial_args = list(name = 'GP_Matern', x_coord_col = 'xname', y_coord_col = 'yname'"
                      )
          stop(msg)
        }
        message("Creating k-means clusters for spatial random effects.")
        if (!is.null(spatial_args[["split_neighbors_by_colname"]])) {
          msg <- paste0("split_neighbors_by_colname argument must be a column in metadata. You may need to specify NULL, or column name in metadata that indicates the 'sample ID'."
                        , " ", spatial_args[["split_neighbors_by_colname"]], " column not found.")
          if (!spatial_args[["split_neighbors_by_colname"]] %in% colnames(metainfo)) {
            stop(msg)
          }
        }
        metainfo <-
        xy_kmeans_clusters(metainfo
                           , x_coord_col = spatial_args[["x_coord_col"]]
                           , y_coord_col = spatial_args[["y_coord_col"]]
                           , k = spatial_args[["k"]]
                           , k_prop_n = spatial_args[["k_prop_n"]]
                           , split_neighbors_by_colname = spatial_args[["split_neighbors_by_colname"]]
                           )
        if (!is.null(spatial_args[["split_neighbors_by_colname"]])) {
          refrmla <- as.formula(paste0("~Matern(1 | "
                                      , spatial_args[["x_coord_col"]]
                                      , "_cluster"
                                      , " + "
                                      , spatial_args[["y_coord_col"]]
                                      , "_cluster"
                                      , " %in% "
                                      , spatial_args[["split_neighbors_by_colname"]]
                                      , ")"
                                      ))
        } else {
          refrmla <- as.formula(paste0("~Matern(1 | "
                                      , spatial_args[["x_coord_col"]]
                                      , "_cluster"
                                      , " + "
                                      , spatial_args[["y_coord_col"]]
                                      , "_cluster"
                                      , ")"
                                      ))
        }
      }
    }

  }

  stopifnot(!is.null(rownames(metadata)) & (cellid_colname %in% names(metainfo)))
  if (!(cellid_colname) %in% names(metainfo)) {
    cellid_colname <- "cell_ID"
    metainfo[[cellid_colname]] <- rownames(metadata)
  }
  if (cellid_colname != "cell_ID" & "cell_ID" %in% names(metainfo)) metainfo[["cell_ID"]] <- NULL
  setnames(metainfo, old = cellid_colname, new = "cell_ID")

  if (is.null(targets)) targets <- rownames(assay_matrix)
  mTerms <- all.vars(formula)
  inla_term <- NULL
  if (!is.null(spatial_model) && spatial_model[["name"]] == "GP_INLA") {
     term_var <- attr(terms(terms(formula)), "variables")
     inla_term <-
     lapply(term_var, function(x) {
       if (is.call(x)) {
         if (any(grepl("INLA|f\\(", as.character(x[[1]])))) {
           unlist(lapply(x[2:length(x)], as.character))
         }
       } else {
        NULL
       }
     })
     inla_term <- unlist(inla_term)
     mTerms <- setdiff(mTerms, inla_term)
  }
  if (!is.null(spatial_model) && spatial_model[["name"]] == "GP_Matern") {
     spamm_terms <- all.vars(spatial_args[["spatial_random_effect"]])
     xycols <- c(spatial_args[["x_coord_col"]], spatial_args[["y_coord_col"]])
     mTerms <- c(mTerms, spamm_terms)
     mTerms <- c(mTerms, xycols)
  }
  mTerms <- setdiff(mTerms, "1")

  if (!is.null(groupVar)) {
    # check if groupVar is in model formula terms
    if (!groupVar %in% mTerms) {
      stop("Error: groupVar needs to be defined as fixed effect in the model.\n")
    }
  }

  Wmat <- NULL
  celltype_ref <- NULL
  if (!missing(neighborhood_counts)) {
    warning("'neighborhood_counts' argument is deprecated and will be removed in the future. Use `pre_de_obj` argument instead.")
    stopifnot("neighborhood_counts must be NULL or an object of neighborexpr class returned by 'measure_neighbor_expr_by_celltype' function" = inherits(neighborhood_counts, "neighborexpr"))

    if (length(setdiff(colnames(assay_matrix), metainfo[[cellid_colname]])) > 0 ||
       length(setdiff(metainfo[[cellid_colname]], colnames(assay_matrix))) > 0
    ) {
      warning("Not identical set of cells between assay matrix and meta.data")
      mcells <- nrow(metainfo)
      acells <- ncol(assay_matrix)
      comm <- length(intersect(colnames(assay_matrix), metainfo[[cellid_colname]]))
      message(paste0(mcells, " cells in metadata."))
      message(paste0(acells, " cells in assay_matrix."))
      message(paste0(comm, " cells common between metadata and assay."))
    }

  } else if (!is.null(pre_de_obj)) {
    stopifnot("pre_de_obj must be NULL or an object of prede class returned by 'pre_de' function" = inherits(pre_de_obj, "prede"))
    neighborhood_counts <- pre_de_obj$nblist
    if (missing(neighbor_expr_cell_type_metadata_colname)) {
      stop("When using `pre_de_obj`, must provide a `neighbor_expr_cell_type_metadata_colname` argument for calculating neighbor celltype expression.")
    }
    stopifnot("provided `neighbor_expr_cell_type_metadata_colname` not found in metadata." =
                neighbor_expr_cell_type_metadata_colname %in% names(metainfo))
    neighbor_expr_overlap_agg <- match.arg(neighbor_expr_overlap_agg)
    if (!is.null(neighbor_expr_overlap_weight_colname)) {
      stopifnot(neighbor_expr_overlap_weight_colname %in% colnames(pre_de_obj$cell_adjacency_dt))
    }
    celltype_ref <- pre_de_obj$cell_adjacency_dt[from == to][, c("from", paste0(neighbor_expr_cell_type_metadata_colname, "_from")), with = FALSE]
    setnames(celltype_ref, paste0(neighbor_expr_cell_type_metadata_colname, "_from"), "ctcol__")

    Wmat <- make_W(data.table(cell_ID = colnames(assay_matrix)) ## manifest of all cells, even if not adjacent to any others
                   , pre_de_obj$cell_adjacency_dt[from != to]
                   , neighbor_expr_cell_type_metadata_colname
                   , neighbor_expr_overlap_agg
                   , "cell_ID"
                   , neighbor_expr_overlap_weight_colname
    )
    ### need a check here?
    Wmat <- Wmat[colnames(assay_matrix), colnames(assay_matrix)]

    if (neighbor_expr_totalcount_normalize) {
      if (is.null(neighbor_expr_totalcount_scalefactor)) {
        message("calculating scalefactors using assay_matrix for totalcount-normalizing neighbor_expr covariates")
        neighbor_expr_totalcount_scalefactor <- Matrix::colSums(assay_matrix)
        mean_scalefactor <- mean(neighbor_expr_totalcount_scalefactor)
        neighbor_expr_totalcount_scalefactor[neighbor_expr_totalcount_scalefactor == 0] <- 1
        neighbor_expr_totalcount_scalefactor <- mean_scalefactor / neighbor_expr_totalcount_scalefactor
      } else {
        stopifnot("provided neighbor_expr_totalcount_scalefactor length does not match number of cells in assay_matrix." =
                    length(neighbor_expr_totalcount_scalefactor) == ncol(assay_matrix))
        stopifnot("provided neighbor_expr_totalcount_scalefactor should be a numeric vector of elements corresponding to each cell in `assay_matrix`." =
                    is.numeric(neighbor_expr_totalcount_scalefactor))
        stopifnot("provided neighbor_expr_totalcount_scalefactor should be a named numeric vector with names corresponding to the cell IDs in `assay_matrix`." =
                   !is.null(names(neighbor_expr_totalcount_scalefactor)))
        stopifnot("provided neighbor_expr_totalcount_scalefactor should be a named numeric vector with names corresponding to the cell IDs in `assay_matrix`." =
                   all(colnames(assay_matrix) %in% names(neighbor_expr_totalcount_scalefactor)))
      }
    }

  }
  # check if terms in model are in sData
  if (is.null(neighborhood_counts)) neighborhood_counts <- list()
  neighbor_terms <- setdiff(mTerms, names(metainfo))
  extract_names <- unique(gsub("_expr$|_intensity$|_cellcount$", "", neighbor_terms))
  candidate_neighbor_terms <- c()
  if (length(neighborhood_counts) > 0) {
    ## potentially lighten the list of matrices passed through to deFunc
    neighborhood_counts$neighbor_expr_byct <- neighborhood_counts$neighbor_expr_byct[extract_names]
    candidate_neighbor_terms <- unlist(lapply(names(neighborhood_counts$neighbor_expr_byct)
                                      , function(x) paste0(x, c("", "_intensity", "_expr", "_cellcount"))))

  } else if (!is.null(pre_de_obj)) {
    candidate_neighbor_terms <- unlist(lapply(c("allct", "otherct", unique(celltype_ref[["ctcol__"]]))
                                      , function(x) paste0(x, c("", "_intensity", "_expr", "_cellcount"))))
    candidate_neighbor_terms <- gsub("\\-|\\ ", "_", candidate_neighbor_terms)
    if (length(setdiff(colnames(assay_matrix), metainfo[[cellid_colname]])) == 0 &
       length(candidate_neighbor_terms > 0)) {
      warning(paste0("The cells in provided `assay_matrix` are the same as those in the `metadata`, "
                     , "and a `pre_de_obj` was passed.\n"
                     , "typical usage would be: smi_de(assay_matrix = FULL_assay_matrix, metadata = CELLSTOANALYZE_metadata)"
                     , "Warning that the `assay_matrix` used to calculate gene expression in neighboring cells should typically contain all of the neighbors of the cells in the metadata.\n"
                     )
             )
    }
  }
  neighbor_terms <- intersect(neighbor_terms, candidate_neighbor_terms)
  mTerms <- setdiff(mTerms, neighbor_terms)

  missingTerms <- setdiff(c(mTerms, neighbor_terms)
                          , c(names(metainfo)
                             , names(neighborhood_counts$neighbor_expr_byct)
                             , candidate_neighbor_terms)
                         )
  if (length(missingTerms) > 0) {
    stop(paste0("Error: ", paste0(missingTerms, collapse = ", "), "were not found in in the meta.data slot of the passed seurat object.\n"))
  }
  pDat <- metainfo[, c("cell_ID", mTerms), with = FALSE]

  ### Add check to ensure that covariates in model have more than one unique value.
  uniq_check <- unlist(pDat[, lapply(.SD, uniqueN), .SDcols = c(mTerms)])
  non_uniq <- names(uniq_check)[which(uniq_check == 1)]
  if (length(non_uniq) > 0) {
    warning(paste0("covariates: ", paste0(non_uniq, collapse = ", "), " in formula have only one unique value in the dataset.\n"))
    for (ii in seq_along(non_uniq)) {
      warning(paste0(non_uniq[ii], " unique value: ", pDat[1][[non_uniq[ii]]]))
    }
    message(paste0("(", nrow(pDat), " total cells)"))
    stop("Non-unique covariates will cause error in regression models, such as:\n"
         , "'contrasts can only be applied to factors with two or more levels'")
  }
  if (!is.null(groupVar)) {
    if (!is.numeric(pDat[[groupVar]])) {
      if (missing(groupVar_levels)) {
        pDat[[groupVar]] <- as.factor(pDat[[groupVar]])
      } else {
        pDat[[groupVar]] <- factor(pDat[[groupVar]], levels = groupVar_levels)
      }
    }
  }

  for (i in setdiff(names(pDat), c("cell_ID", groupVar))) {
    if (inherits(pDat[[i]], "character")) {
      pDat[, i] <- as.factor(pDat[[i]])
    }
  }
  updatedFormula <- formula(paste("y", as.character(formula)[2]
                                  , sep = " ~ "))

  if (!isTRUE(verbose)) {
    message("Fitting model to targets")
  }

  if (nCores > 1) {
    assay_expr <- new.env()
    assay_expr$assay_expr <- assay_matrix
    if (Sys.info()["sysname"] != "Windows") {
      mixedOut <- parallel::mclapply(targets
                                     , deFunc
                                     , groupVar
                                     , groupVar_levels
                                     , pDat
                                     , updatedFormula
                                     , family = family
                                     , assay_expr
                                     , neighbor_expr_overlap_weight_colname
                                     , neighbor_expr_overlap_agg
                                     , neighbor_expr_cell_type_metadata_colname
                                     , neighbor_expr_totalcount_normalize
                                     , neighbor_expr_totalcount_scalefactor
                                     , Wmat
                                     , celltype_ref
                                     , neighborhood_counts
                                     , neighbor_terms
                                     , typ = "parallel"
                                     , spatial_args
                                     , verbose = verbose
                                     , mc.cores = nCores
                                     , ...)
    } else {
      cl <- parallel::makeCluster(getOption("cl.cores", nCores))
      on.exit(try(parallel::stopCluster(cl), silent = TRUE), add = TRUE) # NEW!!!
      mixedOut <- parallel::parLapply(cl
                                      , targets
                                      , deFunc
                                      , groupVar
                                      , pDat
                                      , updatedFormula
                                      , family
                                      , assay_expr
                                      , neighbor_expr_overlap_weight_colname
                                      , neighbor_expr_overlap_agg
                                      , neighbor_expr_cell_type_metadata_colname
                                      , neighbor_expr_totalcount_normalize
                                      , neighbor_expr_totalcount_scalefactor
                                      , Wmat
                                      , celltype_ref
                                      , neighborhood_counts
                                      , neighbor_terms
                                      , typ
                                      , spatial_args
                                      , verbose = verbose
                                      , ...)
      suppressWarnings(parallel::stopCluster(cl)) # redundant with on.exit, but can be left in to ensure cluster is stopped in case of error in parLapply
    }
  } else {
    mixedOut <- lapply(targets
                       , deFunc
                       , groupVar
                       , groupVar_levels
                       , pDat
                       , updatedFormula
                       , family = family
                       , assay_matrix
                       , neighbor_expr_overlap_weight_colname
                       , neighbor_expr_overlap_agg
                       , neighbor_expr_cell_type_metadata_colname
                       , neighbor_expr_totalcount_normalize
                       , neighbor_expr_totalcount_scalefactor
                       , Wmat
                       , celltype_ref
                       , neighborhood_counts
                       , neighbor_terms
                       , typ = "non-parallel"
                       , spatial_args
                       , verbose = verbose
                       , ...)
  }
  names(mixedOut) <- targets

  return_obj <- list(results = mixedOut
                     , groupVar = groupVar
                     , modelterms = mixedOut[[1]]$terms
                     , targets = targets
                     )

  class(return_obj) <- append(class(return_obj), "smide")
  return(return_obj)
}

deFunc <- function(target, groupVar, groupVar_levels, pDat
                   , formula, family
                   , assay_expr
                   , neighbor_expr_overlap_weight_colname
                   , neighbor_expr_overlap_agg
                   , neighbor_expr_cell_type_metadata_colname
                   , neighbor_expr_totalcount_normalize = TRUE
                   , neighbor_expr_totalcount_scalefactor = NULL
                   , Wmat
                   , celltype_ref
                   , neighbor_counts
                   , neighbor_terms
                   , typ
                   , spatial_args
                   , verbose = TRUE
                   , ...) {
  if (isTRUE(verbose)) {
    message_parallel(paste0("Fitting model to target: ", target))
  }
  if (length(neighbor_terms) > 0) {
    extract_names <- unique(gsub("_expr$|_intensity$|_cellcount$", "", neighbor_terms))
    if (length(neighbor_counts) > 0) {
      ### neighbor expression was pre-calculated
      neighb_dt <- extract_neighborexpr_by_gene(neighbor_counts$neighbor_expr_byct
                                                , target
                                                , expr_list_names = extract_names
                                                , cell_IDs = pDat[["cell_ID"]])

      pDat <- cbind(pDat, neighb_dt)
      neighbor_ct_count <- data.table::copy(neighbor_counts$adjacency_counts_by_ct)
      names(neighbor_ct_count) <- gsub("\\-|\\ ", "_", names(neighbor_ct_count))

      pDat <- merge(pDat, neighbor_ct_count, suffixes = c("_expr", "_cellcount"))
      all_celltypes <- setdiff(names(neighbor_ct_count), "cell_ID")
      ref_celltype <- neighbor_counts$ref_celltype
      other_ctypes <- setdiff(all_celltypes
                              , c(ref_celltype, "allct", "otherct"))
      # pDat[, otherct_cellcount:=rowSums(.SD), .SDcols = paste0(other_ctypes, "_cellcount")]
      # pDat[, allct_cellcount:=rowSums(.SD), .SDcols = paste0(all_celltypes, "_cellcount")]
      for (cc in intersect(c("allct", "otherct"), names(pDat))) {
        setnames(pDat, old = cc, new = paste0(cc, "_expr"))
      }
    } else {
      ### calculate the neighbor expression on the fly
      if (typ == "parallel") {
        nex <- assay_expr$assay_expr[target, rownames(Wmat)]
      } else {
        nex <- assay_expr[target, rownames(Wmat)]
      }
      names(nex) <- rownames(Wmat)
      if (neighbor_expr_totalcount_normalize) {
        nex <- nex * neighbor_expr_totalcount_scalefactor[rownames(Wmat)]
      }
      # rm(assay_matrix); gc()
      extract_names <- unique(gsub("_expr$|_intensity$|_cellcount$", "", neighbor_terms))
      uniq_celltypes <- unique(celltype_ref[["ctcol__"]])
      ind_ct_to_calc <- setdiff(extract_names, c("otherct", "allct"))

      modeled_celltypes <- celltype_ref[match(pDat[["cell_ID"]], from)]
      ref_celltype <- unique(modeled_celltypes[["ctcol__"]])

      nblist <- vector(mode = "list", length = length(ref_celltype))
      ict <- 1L
      for (ict in seq_along(ref_celltype)) {
          neighbor_expr_byct <- vector(mode = "list", length = length(extract_names))
          names(neighbor_expr_byct) <- gsub("\\-|\\ ", "_", extract_names)
          ct_ref <- celltype_ref[ctcol__ == ref_celltype[ict], from]
          if ("allct" %in% extract_names) {
            neighbor_expr_byct[["allct"]] <-
              Wmat[ct_ref, ] %*% nex
          }
          if ("otherct" %in% extract_names) {
            ctidx <- setdiff(colnames(Wmat), ct_ref)
            neighbor_expr_byct[["otherct"]] <-
              Wmat[ct_ref, ctidx] %*% nex[ctidx]
          }
          if (length(ind_ct_to_calc) > 0) {
            for (xx in ind_ct_to_calc) {
              ctidx <- celltype_ref[ctcol__ == xx, from]
              neighbor_expr_byct[[gsub("\\-|\\ ", "_", xx)]] <-
                Wmat[ct_ref, ctidx, drop = FALSE] %*%
                nex[ctidx]
            }
          }
          nblist[[ict]] <- do.call(cbind, neighbor_expr_byct)
          #### Work in neighbor cellcount calculations in here.. TO-DO
          #####
      }
      neighb_dt <- as.data.table(do.call(rbind, nblist)[pDat[["cell_ID"]], ])
      colnames(neighb_dt) <- paste0(names(neighbor_expr_byct), "_expr")
      pDat <- cbind(pDat, neighb_dt)
    }
  }

  if (is.character(family)) {
    if (family == "beta") {
      family <- "beta_family"
      warning("please use ", sQuote("beta_family()"), " rather than ",
              sQuote("\"beta\""), " to specify a Beta-distributed response")
    }
    family_char <- family
  } else {
    family_char <- family()$family
  }

  if (typ == "parallel") {
    y <- matrix(assay_expr$assay_expr[target, ], ncol = 1)
    rownames(y) <- colnames(assay_expr$assay_expr)
  } else {
    y <- matrix(assay_expr[target, ], ncol = 1)
    rownames(y) <- colnames(assay_expr)
  }
  if (grepl("(nbinom|pois)", family_char)) {
    if (any(abs(y[, 1] - round(y[, 1])) > 0.001)) {
      warning(sprintf("non-integer counts in a %s model",
                      family_char))
      message_parallel(paste0("converting counts to nearest integer"))
      y[, 1] <- round(y[, 1])
    }
  }
  colnames(y) <- "y"
  dat <- merge(data.table(cell_ID = rownames(y), y = y[, 1])
               , pDat, by = c("cell_ID"))

  re_terms <- lme4::findbars(formula)
  has_re <- length(re_terms) > 0 || !is.null(spatial_args)
  stopifnot("Only 1 random effect grouping level currently supported" =
              length(re_terms) <= 1)
  convergence_error <- FALSE
  err <- NULL
  model_warning_msg <- emmeans_warning_msg <- NULL
  fixed_formula <- NULL ## only needed in nebula::nebula case; remains null for fixed effect only models
  mod_time <- system.time({
  if (!has_re) { ## if no random effect, either glm or MASS::glm.nb (poisson, nb, or gaussian families)
    if (family_char == "nbinom2") {
      fittype <- "MASS::glm.nb"
      mod <-
      withCallingHandlers({
        tryCatch({
          MASS::glm.nb(formula
                              , data = dat
                              , ...
                              )
        }, error = function(e) {
            convergence_error <<- TRUE
            err <<- conditionMessage(e)
            "convergence_error"
          })
      }, warning = function(w) {
        model_warning_msg <<- conditionMessage(w)
        invokeRestart("muffleWarning")
      })

      if (!convergence_error) model_warning_msg <- paste0("converged = ", mod$converged, " ", model_warning_msg)
    } else {
      fittype <- "stats::glm"
      mod <-
      withCallingHandlers({
        tryCatch({
          stats::glm(formula
                    , data = dat
                    , family = family_char
                    , ...
                    )
          }, error = function(e) {
            convergence_error <<- TRUE
            err <<- conditionMessage(e)
            "convergence_error"
          })
      }, warning = function(w) {
        model_warning_msg <<- conditionMessage(w)
        invokeRestart("muffleWarning")
      })
    }
  } else {
    re_formula <- NULL
    if (is.null(spatial_args)) {
      re_formula <- as.formula(paste0("~", paste0(re_terms, collapse = "+")))
    }
    tt <- terms(formula)
    re_vars <- lapply(re_terms, function(x) gsub("1 | ", "", x))
    re_vars <- unlist(lapply(seq_along(re_terms), function(ii) gsub("^.*\\|[\\ ]+", "", re_terms[ii])))
    re_frmlas <- unlist(lapply(seq_along(re_terms), function(ii) gsub("\\|.*$", "", re_terms[ii])))
    re_frmlas <- lapply(re_frmlas, function(x) {
      as.formula(paste0("~", x))
    })
    names(re_frmlas) <- re_vars
    allvar <- as.character(attr(tt, "variables"))[-1]
    response <- allvar[attr(tt, "response")]
    offsetv <- allvar[attr(tt, "offset")]
    allterms <- attr(terms.formula(formula), "term.labels")
    fixed_terms <- setdiff(allterms, re_terms)
    if (length(fixed_terms) == 0) fixed_terms <- "1"
    if (length(offsetv) > 0) fixed_terms <- c(fixed_terms, paste0("offset(", offsetv, ")"))
    fixed_formula <-  as.formula(paste0(response, " ~ ", paste0(fixed_terms, collapse = "+")))

    ### spatial random effect models
    if (!is.null(spatial_args)) {
       if (spatial_args[["name"]] == "GP_Matern") {
         fittype <- "spaMM::fitme"
         newfrmla <- update.formula(formula, paste0(".~.+", labels(terms(spatial_args[["spatial_random_effect"]]))))
         non_model_args <- c("spatial_random_effect", "k", "k_prop_n", "x_coord_col", "y_coord_col", "split_neighbors_by_colname")
         extra_args <- spatial_args[c(setdiff(names(spatial_args), non_model_args))]

         famchar <- switch(family
                           , gaussian = "gaussian"
                           , nbinom2 = "negbin"
                           , poisson = "poisson")
         mod <-
         withCallingHandlers({
           tryCatch({
             capture_output(
             do.call(spaMM::fitme
                     , c(list(formula = newfrmla, data = dat, family = famchar), extra_args))

             )$result

             }, error = function(e) {
               convergence_error <<- TRUE
                err <<- conditionMessage(e)
               "convergence_error"
           })
          }, warning = function(w) {
            model_warning_msg <<- conditionMessage(w)
            invokeRestart("muffleWarning")
          })
       }
       if (spatial_args[["name"]] == "GP_INLA") {
         fittype <- "INLA::inla"
         modelmat <- model.matrix(formula, data = dat)
         realnames <- colnames(modelmat)
         modeldt <- as.data.table(modelmat)
         inlanames <- paste0("Xf", 1:ncol(modeldt))
         data.table::setnames(modeldt, names(modeldt), inlanames)
         newfrmla <- as.formula(paste0("y ~ -1+ ", paste0(colnames(modeldt), collapse = "+")))
         newfrmla <- update.formula(newfrmla, paste0(".~.+", labels(terms(spatial_args[["spatial_random_effect"]]))))

         offsetcol <- NULL
         offsetuse <- NULL
         if (length(offsetv) > 0) {
           offsetuse <- eval(parse(text = offsetv), envir = as.data.frame(dat))
         } else {
           offsetuse <- NA
         }

         plist <- list(s = 1:ncol(spatial_args[["A"]]))
         for (nm in names(modeldt)) {
           plist[[nm]] <- modeldt[[nm]]
         }
         inlastack <- INLA::inla.stack(
           tag = "est"
           , data = list(dat[, .(y)]) #should include outcome and offset, but not covariates
           , effects = plist
           , a = c(spatial_args[["A"]]
                  , lapply(1:(length(plist) - 1), function(x) 1)
                  )
         )
         lclists <- prepare_lincombs(dat = dat, fitmethod = "INLA::inla", original_formula = formula, fam = family_char, inlanames = inlanames)
         famchar <- switch(family_char
                           , gaussian = "gaussian"
                           , nbinom2 = "nbinomial"
                           , poisson = "poisson"
                           )

         non_model_args <- c("x_coord_col"
                              , "y_coord_col"
                              , "max.edge.mesh"
                              , "cutoff.mesh"
                              , "offset.mesh"
                              , "alpha.inlaprior"
                              , "prior.range.inlaprior"
                              , "prior.sigma.inlaprior"
                              , "constr.inlaprior"
                              , "mesh"
                              , "spatial_random_effect", "priors", "name", "A")
         extra_args <- spatial_args[c(setdiff(names(spatial_args), non_model_args))]
         mod <-
         withCallingHandlers({
           tryCatch({
         do.call(INLA::inla
                   , c(list(formula = newfrmla
                     , e = offsetuse
                     , data = INLA::inla.stack.data(inlastack, spde = spatial_args[["priors"]])
                        , control.predictor = list(A = INLA::inla.stack.A(inlastack)
                                                 , compute = TRUE)
                        , family = famchar
                        , lincomb = c(lclists[["emm_lc_list"]]
                                    , lclists[["pairwise_lc_list"]]
                                    , lclists[["onevrest_lc_list"]]
                                    , lclists[["onevall_lc_list"]])
                        , control.inla = list(int.strategy = "eb")
                 ), extra_args)
                 )
           }, error = function(e) {
             convergence_error <<- TRUE
              err <<- conditionMessage(e)
             "convergence_error"
         })
        }, warning = function(w) {
          model_warning_msg <<- conditionMessage(w)
          invokeRestart("muffleWarning")
        })

         mod$family <- family
       }
    } else {
      if (family_char == "gaussian") {
        fittype <- "lme4::lmer"
        mod <-
        withCallingHandlers({
         tryCatch({
           lme4::lmer(formula
                      , data = dat
                      , ...
                      )
           }, error = function(e) {
             convergence_error <<- TRUE
              err <<- conditionMessage(e)
             "convergence_error"
         })
        }, warning = function(w) {
          model_warning_msg <<- conditionMessage(w)
          invokeRestart("muffleWarning")
        })
      } else if (family_char %in% c("nbinom2", "poisson")) {
        fittype <- "nebula::nebula"
        modelv <- "NBGMM" # neg binomial gamma mixed model
        if (family_char == "poisson") modelv <- "PMM" # poisson gamma mixed model
        data.table::setkeyv(dat, re_vars)
        offsetcol <- NULL
        offsetuse <- NULL
        if (length(offsetv) > 0) {
          offsetuse <- eval(parse(text = offsetv), envir = as.data.frame(dat))
        }
        mod <-
        withCallingHandlers({
          tryCatch({
            #suppressMessages(
             capture_output(
              nebula::nebula(
                count = as(matrix(dat[["y"]], nrow = 1), "dgCMatrix")
                , id = dat[[re_vars]]
                , pred = model.matrix(fixed_formula, data = dat)
                , offset = offsetuse
                , covariance = TRUE
                , ncore = 1
                , model = modelv
              )
            )$result
            ## attach overdispersions if needed
          }, error = function(e) {
            convergence_error <<- TRUE
            err <<- conditionMessage(e)
            "convergence_error"
          })
        }, warning = function(w) {
          model_warning_msg <<- conditionMessage(w)
          invokeRestart("muffleWarning")
        })

      }  else {
        fittype <- "MASS::glmmPQL"
        mod <-
        withCallingHandlers({
          tryCatch({
            MASS::glmmPQL(fixed = fixed_formula
                          , random = re_formula
                          , family = family
                          , data = dat
                          , ...
                          )
            }, error = function(e) {
              convergence_error <<- TRUE
              err <<- conditionMessage(e)
              "convergence_error"
         })
         }, warning = function(w) {
           model_warning_msg <<- conditionMessage(w)
           invokeRestart("muffleWarning")
         })
      }
    }
  }
  })

  if (fittype != "INLA::inla") {
    mod_out <-
      withCallingHandlers({
                  summarize_model(mod
                                 , modtime = mod_time
                                 , groupVar
                                 , dat
                                 , fittype
                                 , original_formula = formula
                                 , fixed_formula = fixed_formula
                                 , error_msg = err
                                 )
         }, warning = function(w) {
           emmeans_warning_msg <<- conditionMessage(w)
           invokeRestart("muffleWarning")
         })

  }
  if (!is.null(spatial_args)) {
    if (fittype == "INLA::inla") {
      mod_out <-
      withCallingHandlers({
        postprocess_inla_contrasts(lclists, mod, inlanames, realnames, mod_time, error_msg = err)
      }, warning = function(w) {
        emmeans_warning_msg <<- conditionMessage(w)
        invokeRestart("muffleWarning")
      })

      re_name <- all.vars(spatial_args[["spatial_random_effect"]])[[1]]
      mod_out[["spatial_random_effect"]] <- cbind(data.table(x = spatial_args[["mesh"]]$loc[, 1], y = spatial_args[["mesh"]]$loc[, 2])
                                                  , mod$summary.random[[re_name]])
      mod_out[["spatial_random_effect"]][, target := target]

    }
    if (fittype == "spaMM::fitme") {
      joiners <- all.vars(spatial_args[["spatial_random_effect"]])
      xycols <- c(spatial_args[["x_coord_col"]], spatial_args[["y_coord_col"]])
      mod_out[["spatial_random_effect"]] <-
      cbind(data.table(re = spaMM::ranef(mod)[[1]], cluster = names(spaMM::ranef(mod)[[1]]))
        , dat[, head(.SD, 1), by = c(joiners)][, joiners, with = FALSE]
      )
      mod_out[["spatial_random_effect"]] <- merge(dat[, c(xycols, joiners), with = FALSE], mod_out[["spatial_random_effect"]]
                      , by = joiners
                      , sort = FALSE)
      mod_out[["spatial_random_effect"]][, target := target]
    }
    if (fittype == "nebula::nebula") {

    }
  }

  mod_out[["model_summary"]][, target := target]
  for (dt in mod_out[["pairwise_list"]]) dt[, target := target]
  for (dt in mod_out[["onevall_list"]]) dt[, target := target]
  for (dt in mod_out[["onevrest_list"]]) dt[, target := target]
  for (dt in mod_out[["emm_list"]]) dt[, target := target]

  if (!is.null(model_warning_msg)) {
    mod_out[["model_summary"]][, msg := model_warning_msg]
  }
  if (!is.null(emmeans_warning_msg)) {
    for (dt in mod_out[["pairwise_list"]]) dt[, msg := emmeans_warning_msg]
    for (dt in mod_out[["onevall_list"]]) dt[, msg := emmeans_warning_msg]
    for (dt in mod_out[["onevrest_list"]]) dt[, msg := emmeans_warning_msg]
    for (dt in mod_out[["emm_list"]]) dt[, msg := emmeans_warning_msg]
  }
  return(mod_out)
}
