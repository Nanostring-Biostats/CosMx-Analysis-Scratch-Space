#' Extract results from DE object by contrast comparison, model variable, gene
#'
#' @param smide_results object of class "smide" returned by call to run_de or smi_de functions.
#' @param comparisons one or more of c("pairwise", "one.vs.rest", "one.vs.all", "emmeans", "model_summary").
#'                    Corresponding results from emmeans package are returned for variables which were included
#'                    as fixed effects in the DE model formula.
#'                    For continuous variables, only the 'emmeans' returns a result,
#'                    corresponding to that of the mean of the untransformed covariate and the corresponding model estimated marginal response of the outcome at that mean.
#' * pairwise : pairwise contrasts between levels of a categorical fixed effect variable in DE model. See pairwise.emmc in ?emmeans::`contrast-methods` .
#' * one.vs.rest : one category vs. average of all other category contrasts, for categorical fixed effect variables.
#'
#'                 This contrast will be empty for categorical variables with only 2 levels,
#'                 as the pairwise contrast is equivalent. The weights argument 'wts' is set to correspond to # of cells in a category,
#'                 so that these contrasts don't emphasize niche or categorical groups with very few cells. See del.eff.emmc in ?emmeans::`contrast-methods`.
#'                 (also see here for more discussion on using weights : \href{https://github.com/rvlenth/emmeans/issues/346}{emmeans github issue}).
#' * one.vs.all : one category vs. average of all categories (including the one being compared to), for categorical fixed effect variables.
#'                 See eff.emmc in ?emmeans::`contrast-methods`.
#'                 This contrast will be empty for categorical variables with only 2 levels. The weights argument 'wts' is set to correspond to # of cells in a category,
#'                 so that these contrasts don't emphasize niche or categorical groups with very few cells.
#' * emmeans : Model estimated marginal means of the response level, at different levels of the predictor.
#'              Ratios of these can be used to estimate fold change. Lots of vignettes and examples available online.
#'              See here, for example (\href{https://cran.r-project.org/web/packages/emmeans/vignettes/}{emmeans-vignettes})
#' * model_summary : summary tables for fixed effects and std deviations of variance components.
#'                  (i.e., this is equal to sqrt(mse) in linear regression models, or the dispersion parameter in a negative binomial model.)
#' * all : if 'all' is contained anywhere in comparisons argument, a list of all available results are returned.
#' @param variable vector of one or more of the fixed effect variables in DE model.  Can be "all".
#' @param targets vector of one or more of the targets analyzed in smide_results object. Can be "all".
#' @return list of data.tables corresponding to targets, comparisons, and fixed effect variables from DE models.
#'
#' @examples
#'
#' library(Seurat)
#' library(data.table)
#' setDTthreads(1)
#' datadir <- system.file("extdata", package = "smiDE")
#' sem <- readRDS(paste0(datadir, "/small_nsclc.rds"))
#' metainfo <- data.table::data.table(sem@meta.data)
#'
#' totalcount_scalefactors <- mean(metainfo[["totalcounts"]]) / metainfo[["totalcounts"]]
#' names(totalcount_scalefactors) <- colnames(sem)
#' sem <- Seurat::SetAssayData(
#'   sem,
#'   "data",
#'   sem[["RNA"]]@counts %*% Matrix::Diagonal(x = totalcount_scalefactors, names = colnames(sem))
#' )
#' pre_de_obj <-
#'   pre_de(
#'     metadata = metainfo,
#'     cell_type_metadata_colname = "cell_type",
#'     split_neighbors_by_colname = "tissue",
#'     mm_radius = 0.05,
#'     sdimx_colname = "sdimx",
#'     sdimy_colname = "sdimy",
#'     verbose = TRUE
#'   )
#'
#' fibroblast_cells <- metainfo[cell_type == "fibroblast", cell_ID]
#' de_results <-
#'   smi_de(
#'     assay_matrix = sem[["RNA"]]@counts,
#'     metadata = metainfo[cell_ID %in% fibroblast_cells],
#'     formula = ~ RankNorm(otherct_expr) + niche + tissue + offset(log(totalcounts)),
#'     groupVar = "niche",
#'     family = "nbinom2",
#'     targets = rownames(sem)[1:2],
#'     pre_de_obj = pre_de_obj,
#'     neighbor_expr_cell_type_metadata_colname = "cell_type",
#'     neighbor_expr_contamination_weight_colname = NULL,
#'     neighbor_expr_contamination_agg = "sum",
#'     neighbor_expr_totalcount_normalize = TRUE,
#'     neighbor_expr_totalcount_scalefactor = totalcount_scalefactors
#'   )
#'
#' results(de_results, comparisons = "pairwise", variable = "niche")[[1]]
#' results(de_results, comparisons = "model_summary", targets = "CD14")[[1]]
#' results(de_results, comparisons = "one.vs.rest", variable = "tissue")[[1]]
#' results(de_results, comparisons = "emmeans", variable = "otherct_expr")
#' results(de_results, comparisons = "emmeans")
#' results(de_results)
#'
#' @export
results <- function(
  smide_results,
  comparisons = c("pairwise", "one.vs.rest", "one.vs.all", "emmeans", "model_summary", "all", "spatial_random_effect"),
  variable = c(smide_results$groupVar, smide_results$modelterms, "all"),
  targets = c("all")
) {
  if ("all" %in% comparisons) comparisons <- c("pairwise", "one.vs.rest", "one.vs.all", "emmeans", "model_summary", "spatial_random_effect")
  if ("all" %in% variable) variable <- smide_results$modelterms

  stopifnot(inherits(smide_results, "smide"))
  stopifnot(
    "'variable' are not all fixed effects terms in model formula (smide_results$modelterms) " =
      all(variable %in% smide_results$modelterms)
  )

  if (length(targets) == 1 && targets == "all") targets <- smide_results$targets
  stopifnot(all(targets %in% smide_results$targets))


  res_names <- c(
    "pairwise" = "pairwise_list",
    "one.vs.rest" = "onevrest_list",
    "one.vs.all" = "onevall_list",
    "emmeans" = "emm_list",
    "model_summary" = "model_summary",
    "spatial_random_effect" = "spatial_random_effect"
  )

  res_out <- smide_results$results[targets]
  outcomp <- list()
  for (comp in setdiff(comparisons, c("model_summary", "spatial_random_effect"))) {
    outd <- list()
    for (vv in variable) {
      outd[[vv]] <- rbindlist(lapply(lapply(res_out, "[[", res_names[comp]), "[[", vv), use.names = TRUE, fill = TRUE)
    }
    outcomp[[comp]] <- rbindlist(outd, use.names = TRUE, fill = TRUE)
  }
  if ("model_summary" %in% comparisons) {
    outcomp[["model_summary"]] <- rbindlist(lapply(res_out, "[[", "model_summary"), use.names = TRUE, fill = TRUE)
  }

  if ("spatial_random_effect" %in% comparisons) {
    outcomp[["spatial_random_effect"]] <- rbindlist(lapply(res_out, "[[", "spatial_random_effect"), use.names = TRUE, fill = TRUE)
  }

  return(outcomp)
}


#' Function used to summarize models with emmeans based contrasts and a model summary.
#'
#' @param mod a model.
#' @param modtime optional. An object produced by 'system.time', if provided it is added to the model_summary in returned object.
#' @param groupVar optional groupVar.
#' If specified, this forces pairwise contrasts to be conducted by this variable, even if there are 20 or more levels.
#' @return a list of lists containing data.tables of emmeans summaries.  Intended to be primarily accessed through the smiDE::results function.
#'
#'
#' @export
#'
summarize_model <- function(
  mod,
  modtime = NULL,
  groupVar = NULL,
  dat = NULL,
  fitmethod,
  original_formula,
  fixed_formula,
  error_msg
) {
  ## output list for:
  ## estimated marginal means (emm_list)
  ## pairwise contrasts (pw_list)
  ## one vs. all groups contrasts (onevall_list)
  ## one vs. the rest contrasts (onevrest_list)
  emm_list <- list()
  pw_list <- list()
  onevall_list <- list()
  onevrest_list <- list()
  allterms <- c()
  if (is.null(error_msg)) {
    linkused <- switch(fitmethod,
      "MASS::glmmPQL" = mod$family$link,
      "stats::glm" = family(mod)$link,
      "lme4::lmer" = family(mod)$link,
      "MASS::glm.nb" = family(mod)$link,
      "nebula::nebula" = "log",
      "spaMM::fitme" = family(mod)$link
    )
  }

  ## gather terms to be summarized
  re_terms <- NULL
  allterms <- fixedterms <- attr(terms.formula(original_formula), "term.labels")
  interactions <- grep(":", attr(terms.formula(original_formula), "term.labels"), value = TRUE)
  if (length(lme4::findbars(original_formula)) > 0) {
    re_terms <- lme4::findbars(original_formula)
    fixedterms <- setdiff(allterms, re_terms)
  }
  fixedterms <- unlist(lapply(fixedterms, function(x) all.vars(as.formula(paste0("~", x)))))
  fixedterms <- unique(fixedterms)
  allterms <- c(fixedterms, interactions) # unlist(union(allterms, fixedterms))


  .is.num <- function(x) {
    nx <- suppressWarnings(as.numeric(as.character(x)))
    !any(is.na(nx))
  }

  ## Fixed effects model terms contrasts
  if (length(fixedterms) > 0) {
    for (term in fixedterms) {
      continuous_var <- FALSE
      if (!is.null(error_msg)) {
        emm_list[[term]] <- data.table(msg = error_msg, term = term)
        pw_list[[term]] <- data.table(msg = error_msg, term = term)
        onevrest_list[[term]] <- data.table(msg = error_msg, term = term)
        onevall_list[[term]] <- data.table(msg = error_msg, term = term)
      } else {
        if (is.factor(dat[[term]])) {
          moduse <- prepare_mod_emmeans(mod, term, fitmethod, fixed_formula, dat, cov_reduce_list = list())
          emob <- suppressMessages(suppressWarnings(emmeans::emmeans(moduse, term, type = "response", data = dat)))
        } else {
          continuous_var <- TRUE
          cov_reduce_list <- make_cov_reduce_list(fixedterms, dat)
          moduse <- prepare_mod_emmeans(mod, term, fitmethod, fixed_formula, dat, cov_reduce_list)
          ## if variable is continuous, evaluate at the marginal mean, and marginal mean + 1 standard deviation
          emob <- suppressMessages(suppressWarnings(emmeans::emmeans(moduse,
            term,
            type = "response",
            data = dat,
            cov.reduce = cov_reduce_list
          )))
        }
        emobdt <- as.data.table(emob)[, ky := "ky"]
        setnames(emobdt, old = c(term), new = c("level"))
        if (!("response" %in% names(emobdt)) && ("emmean" %in% names(emobdt))) {
          setnames(emobdt, old = "emmean", new = "response")
        }
        if (!("response" %in% names(emobdt)) && ("rate" %in% names(emobdt))) {
          setnames(emobdt, old = "rate", new = "response")
        }
        # setnames(emobdt, old=c(names(emobdt)[2], term), new=c("response","level"))
        emobdt[, term := term]
        levs <- emob@levels[[1]]

        if (length(levs) >= 2) {
          wtdt <- as.data.table(emob@grid)
          setnames(wtdt, old = term, new = "level")
          if (length(emob@model.info$nesting[[term]]) > 0) {
            wts <- wtdt[.wgt. > 0, .wgt.]
            names(wts) <- wtdt[.wgt. > 0, level]
          } else {
            wts <- wtdt[, .wgt.]
            names(wts) <- wtdt[, level]
          }
          setnames(wtdt, ".wgt.", "wts")
          levs <- levs[match(names(wts), levs)]
          levdt <- data.table(category = paste0("c", seq_along(levs)))[, level := levs]
          emobdt <- merge(emobdt, levdt, by = "level")
          emobdt <- merge(emobdt, wtdt, by = c("level", intersect(fixedterms, names(wtdt))))
          if (emobdt[, .is.num(level)]) emobdt[, level := paste0(term, "", as.character(level))]

          ### for continuous variables, no meaningful 'total counts' or 'ncells' associated w compared levels
          ### while for categorical variables, add in the observed counts and number of cells for
          ### corresponding categories in the contrast
          stats_all <- dat[, .(level = "avg.all", counts = sum(y), propnz = mean(y > 0), ncells = .N)]
          if (!continuous_var) {
            stats_rest <- rbindlist(lapply(emobdt[["level"]], function(lev) {
              cpc_rest <- dat[dat[[term]] != lev][, .(term_1 = lev, level = "avg.rest", counts = sum(y), propnz = mean(y > 0), ncells = .N)]
              cpc_lev <- dat[dat[[term]] == lev][, .(term_1 = lev, level = lev, counts = sum(y), propnz = mean(y > 0), ncells = .N)]
              rbindlist(list(cpc_rest, cpc_lev))
            }))
            stats_rest <- merge(stats_rest,
              emobdt[, .(term_1 = level, modelest_cpc = response)],
              by = c("term_1")
            )
          }
          if (length(levs) < 10**3 || (!is.null(groupVar) && term == groupVar)) {
            pw <- as.data.table(suppressWarnings(emmeans::contrast(emob, "pairwise", parens = NULL, adjust = "none", data = dat)))
            # if (continuous_var) pw[["contrast"]] <- paste0(paste0(paste0(c("Mean [", "Mean + 1SD ["), term), "]"), collapse=" / ")
          }

          ### one vs. rest and one vs. all contrasts only created if more than 2 levels of the
          ### categorical variable, otherwise they are redundant w/ a pairwise contrast
          if (length(levs) > 2) {
            onevrest <- as.data.table(suppressWarnings(emmeans::contrast(emob, "del.eff", wts = wts, parens = NULL, adjust = "none", data = dat)))
            onevrest[, term_1 := gsub(" effect", "", contrast)]
            onevrest[, contrast := gsub("effect", "vs. avg.rest", contrast)]
            if (any(grepl("^\\(.*\\)$", onevrest[["term_1"]]))) {
              ## if a leading andtrailing parenthesis is added to the contrast label,
              ## homogenize by adding to all levels.
              ## (This can happen in nested case, despite argument parens=NULL above.. may file an issue w/ emmeans to see if a bug.)
              onevrest[!grepl("^\\(.*\\)", term_1), `:=`(
                contrast = gsub("^", "(", gsub("\\ vs\\.", ") vs.", contrast)),
                term_1 = paste0("(", term_1, ")")
              )]
            }
            emobdt[, level_use_join := level]
            stats_rest[, level_use_join := term_1]
            confounded_covariates <- intersect(fixedterms, names(emobdt))
            if (length(confounded_covariates) > 0) {
              for (ii in seq_along(confounded_covariates)) {
                emobdt[["level_use_join"]] <- paste0(emobdt[["level_use_join"]], " ", emobdt[[confounded_covariates[ii]]])
                add_parenth <- onevrest[grepl("\\(", term_1), term_1]
                emobdt[
                  level_use_join %in% gsub("\\(|\\)", "", add_parenth),
                  level_use_join := paste0("(", level_use_join, ")")
                ]
                stats_rest[level != "avg.rest", level_use_join := paste0(level_use_join, " ", emobdt[[confounded_covariates[ii]]])]
                stats_rest[level == "avg.rest", level_use_join := paste0(level_use_join, " ", emobdt[[confounded_covariates[ii]]])]
                stats_rest[
                  level_use_join %in% gsub("\\(|\\)", "", add_parenth),
                  level_use_join := paste0("(", level_use_join, ")")
                ]
              }
            }

            onevrest <- merge(onevrest, emobdt[, c("level_use_join", "response"), with = FALSE],
              by.y = "level_use_join", by.x = "term_1"
            )

            onevall <- as.data.table(suppressWarnings(emmeans::contrast(emob, "eff", wts = wts, parens = NULL, adjust = "none", data = dat)))
            onevall[, term_1 := gsub(" effect", "", contrast)]
            onevall[, contrast := gsub("effect", "vs. avg.all", contrast)]
            if (any(grepl("^\\(.*\\)$", onevall[["term_1"]]))) {
              ## if a leading andtrailing parenthesis is added to the contrast label,
              ## homogenize by adding to all levels.
              ## (This can happen in nested case, despite argument parens=NULL above.. may file an issue w/ emmeans to see if a bug.)
              # onevall[!grepl("^\\(.*\\)", term_1),` :=`(term_1=paste0("(",term_1, ")")
              onevall[!grepl("^\\(.*\\)", term_1), `:=`(
                contrast = gsub("^", "(", gsub("\\ vs\\.", ") vs.", contrast)),
                term_1 = paste0("(", term_1, ")")
              )]
            }


            emobdt[, level_use_join := level]
            if (length(confounded_covariates) > 0) {
              for (ii in seq_along(confounded_covariates)) {
                emobdt[["level_use_join"]] <- paste0(emobdt[["level_use_join"]], " ", emobdt[[confounded_covariates[ii]]])
                add_parenth <- onevall[grepl("\\(", term_1), term_1]
                emobdt[
                  level_use_join %in% gsub("\\(|\\)", "", add_parenth),
                  level_use_join := paste0("(", level_use_join, ")")
                ]
              }
            }
            onevall <- merge(onevall, emobdt[, c("level_use_join", "response"), with = FALSE],
              by.y = "level_use_join", by.x = "term_1"
            )

            emobdt[, level_use_join := NULL]
            if (linkused == "identity") {
              onevrest[, fold_change := response / (response - estimate)]
              onevall[, fold_change := response / (response - estimate)]
            } else if (linkused == "log") {
              onevrest[, fold_change := ratio]
              onevall[, fold_change := ratio]
              pw[, fold_change := ratio]
            } else {
              onevrest[, fold_change := NA_real_]
              onevall[, fold_change := NA_real_]
              pw[, fold_change := NA_real_]
            }
            ## Add in counts and number of cells metrics
            onevrest <- merge(onevrest,
              stats_rest[level != "avg.rest"][, .(term_1 = level_use_join, counts_1 = counts, propnz_1 = propnz, ncells_1 = ncells)],
              by = c("term_1"),
              sort = FALSE
            )

            onevrest <- merge(onevrest,
              stats_rest[level == "avg.rest"][, .(term_1 = level_use_join, counts_2 = counts, propnz_2 = propnz, ncells_2 = ncells)],
              by = c("term_1"),
              sort = FALSE
            )

            onevrest[, modelest_cpc_1 := response]
            onevrest[, modelest_cpc_2 := modelest_cpc_1 / fold_change]

            onevall <- merge(onevall,
              stats_rest[level != "avg.rest"][, .(term_1 = level_use_join, counts_1 = counts, propnz_1 = propnz, ncells_1 = ncells)],
              by = c("term_1"),
              sort = FALSE
            )

            onevall <- merge(onevall,
              stats_all[, .(term_1 = onevall[["term_1"]], counts_2 = counts, propnz_2 = propnz, ncells_2 = ncells)],
              by = c("term_1"),
              sort = FALSE
            )
            #         ,new="pval")
            onevall[, modelest_cpc_1 := response]
            onevall[, modelest_cpc_2 := modelest_cpc_1 / fold_change]
            rm_names <- c("response", "term_1")
            onevall[, (rm_names) := NULL]
            onevrest[, (rm_names) := NULL]

            onevrest_list[[term]] <- rbindlist(list(
              onevrest_list[[term]],
              onevrest[, ` :=`(term = term)]
            ))

            onevall_list[[term]] <- rbindlist(list(
              onevall_list[[term]],
              onevall[, ` :=`(term = term)]
            ))
          }

          ## Manually calculate fold change ratios for models with identity links (i.e., gaussian)
          if (linkused == "identity") {
            fc_pw <-
              merge(emobdt[, c("level", "response", "ky", "category"), with = FALSE],
                emobdt[, c("level", "response", "ky", "category"), with = FALSE],
                by = "ky", allow.cartesian = TRUE, suffixes = c("_1", "_2")
              )
            fc_pw <- fc_pw[category_1 < category_2]
            fc_pw[, fold_change := response_1 / response_2]
            fc_pw[, contrast := paste0(level_1, " / ", level_2)]
            pw[, level_1 := tstrsplit(contrast, " - ")[[1]]]
            pw[, level_2 := tstrsplit(contrast, " - ")[[2]]]
            pw <- merge(pw, fc_pw[, .(level_1, level_2, fold_change)],
              by = c("level_1", "level_2")
            )
            rm_names <- c("level_1", "level_2")
            pw[, (rm_names) := NULL]
          } else if (linkused == "log") {
            pw[, fold_change := ratio]
          } else {
            pw[, fold_change := NA_real_]
          }

          pw[, term_1 := tstrsplit(contrast, "\\ [\\/|-]\\ ")[[1]]]
          pw[, term_2 := tstrsplit(contrast, "\\ [\\/|-]\\ ")[[2]]]
          pw <- merge(pw, emobdt[, .(term_1 = level, modelest_cpc_1 = response)], by = c("term_1"), all.x = TRUE, sort = FALSE)
          pw <- merge(pw, emobdt[, .(term_2 = level, modelest_cpc_2 = response)], by = c("term_2"), all.x = TRUE, sort = FALSE)
          if (!continuous_var) {
            pw <- merge(pw,
              stats_rest[level != "avg.rest"][, .(term_1,
                counts_1 = counts,
                propnz_1 = propnz,
                ncells_1 = ncells
              )],
              by = c("term_1"),
              all.x = TRUE,
              sort = FALSE
            )
            pw <- merge(pw,
              stats_rest[level != "avg.rest"][, .(
                term_2 = term_1,
                counts_2 = counts,
                propnz_2 = propnz,
                ncells_2 = ncells
              )],
              by = c("term_2"),
              all.x = TRUE,
              sort = FALSE
            )
          }

          pw[, ` :=`(term_1 = NULL, term_2 = NULL)]
          if (continuous_var) {
            pw[, l1o := tstrsplit(contrast, "\\ [\\/-]\\ ")[[1]]]
            pw[, l2o := tstrsplit(contrast, "\\ [\\/-]\\ ")[[2]]]
            pw[, l1 := gsub(term, paste0(term, " "), l1o)]
            pw[, l2 := gsub(term, paste0(term, " "), l2o)]
            pw[, l1 := paste0(term, " ", formatC(as.numeric(tstrsplit(l1, " ")[[2]])))]
            pw[, l2 := paste0(term, " ", formatC(as.numeric(tstrsplit(l2, " ")[[2]])))]
            pw[, contrast := gsub(l1o, l1, contrast)]
            pw[, contrast := gsub(l2o, l2, contrast)]
            pw[, ` :=`(l1 = NULL, l2 = NULL, l1o = NULL, l2o = NULL)]
          }
          pw_list[[term]] <- rbindlist(list(
            pw_list[[term]],
            pw[, ` :=`(term = term)]
          ), use.names = TRUE, fill = TRUE)
        }
        emobdt[, ky := NULL]
        emm_list[[term]] <- rbindlist(list(
          emm_list[[term]],
          emobdt[, `:=`(term = term)]
        ))
      }
    }
  }
  if (length(interactions) > 0) {
    for (term in interactions) {
      for (conditionvar in strsplit(term, ":")[[1]]) {
        emvar <- setdiff(strsplit(term, ":")[[1]], conditionvar)
        if (!is.null(error_msg)) {
          emm_list[[term]] <- data.table(msg = error_msg, term = term)
          pw_list[[term]] <- data.table(msg = error_msg, term = term)
          onevrest_list[[term]] <- data.table(msg = error_msg, term = term)
          onevall_list[[term]] <- data.table(msg = error_msg, term = term)
        } else {
          moduse <- prepare_mod_emmeans(mod, term, fitmethod, fixed_formula, dat)
          emob <- suppressMessages(suppressWarnings(emmeans::emmeans(moduse, emvar, type = "response", data = dat, by = conditionvar)))
          emobdt <- as.data.table(emob)[, ky := "ky"]
          setnames(emobdt, old = c(emvar), new = c("level"))
          if (!("response" %in% names(emobdt)) && ("emmean" %in% names(emobdt))) {
            setnames(emobdt, old = "emmean", new = "response")
          }
          if (!("response" %in% names(emobdt)) && ("rate" %in% names(emobdt))) {
            setnames(emobdt, old = "rate", new = "response")
          }
          # setnames(emobdt, old=c(names(emobdt)[2], term), new=c("response","level"))
          emobdt[, term := term]
          levs <- emob@levels[[1]]
          stats_all <- dat[, .(level = "avg.all", counts = sum(y), propnz = mean(y > 0), ncells = .N), by = conditionvar]
          stats_all <- merge(stats_all[, ky := "ky"], data.table(term_1 = unique(levs), ky = "ky"), allow.cartesian = TRUE)[, ky := NULL]
          stats_rest <- rbindlist(lapply(unique(emobdt[["level"]]), function(lev) {
            cpc_rest <- dat[dat[[emvar]] != lev][, .(term_1 = lev, level = "avg.rest", counts = sum(y), propnz = mean(y > 0), ncells = .N), by = c(conditionvar)]
            cpc_lev <- dat[dat[[emvar]] == lev][, .(term_1 = lev, level = lev, counts = sum(y), propnz = mean(y > 0), ncells = .N), by = c(conditionvar)]
            rbindlist(list(cpc_rest, cpc_lev))
          }))
          if (length(levs) >= 2) {
            wtdt <- as.data.table(emob@grid)
            setnames(wtdt, old = emvar, new = "level")
            wtdt[, idx := seq_len(.N)]
            wtdt_split <- split(wtdt, by = conditionvar)
            for (jj in seq_along(wtdt_split)) {
              if (length(emob@model.info$nesting[[term]]) > 0) {
                wts <- wtdt_split[[jj]][.wgt. > 0, .wgt.]
                names(wts) <- wtdt_split[[jj]][.wgt. > 0, level]
              } else {
                wts <- wtdt_split[[jj]][, .wgt.]
                names(wts) <- wtdt_split[[jj]][, level]
              }
              # wts <- wtdt_split[[jj]][.wgt. > 0,.wgt.]
              # names(wts) <- wtdt_split[[jj]][.wgt. > 0,level]
              setnames(wtdt_split[[jj]], ".wgt.", "wts")
              levs <- levs[match(names(wts), levs)]
              levdt <- data.table(category = paste0("c", seq_along(levs)))[, level := levs]
              emobdt_jj <- merge(emobdt, levdt, by = "level")
              emobdt_jj <- merge(emobdt_jj, wtdt_split[[jj]], by = c("level", intersect(fixedterms, names(wtdt))))
              emobdt_jj <- emobdt_jj[match(levs, level)]
              if (emobdt_jj[, .is.num(level)]) emobdt_jj[, level := paste0(term, " ", as.character(level))]
              if (length(levs) < 10**4 || (!is.null(groupVar) && term == groupVar)) {
                pw <- as.data.table(suppressWarnings(emmeans::contrast(emob[emobdt_jj[, idx]], "pairwise", parens = NULL, adjust = "none", data = dat, by = conditionvar)))
              }
              if (length(levs) > 2) {
                onevrest <- as.data.table(suppressWarnings(emmeans::contrast(emob[emobdt_jj[, idx]], "del.eff", wts = emobdt_jj[, wts], parens = NULL, adjust = "none", data = dat, by = conditionvar)))
                onevrest[, term_1 := gsub(" effect", "", contrast)]
                onevrest[, contrast := gsub("effect", "vs. avg.rest", contrast)]
                emobdt_jj[, level_use_join := level]
                confounded_covariates <- setdiff(intersect(fixedterms, names(emobdt)), conditionvar)
                if (length(confounded_covariates) > 0) {
                  for (ii in seq_along(confounded_covariates)) {
                    emobdt_jj[["level_use_join"]] <- paste0(emobdt_jj[["level_use_join"]], " ", emobdt_jj[[confounded_covariates[ii]]])
                    add_parenth <- onevrest[grepl("\\(", term_1), term_1]
                    emobdt_jj[
                      level_use_join %in% gsub("\\(|\\)", "", add_parenth),
                      level_use_join := paste0("(", level_use_join, ")")
                    ]
                  }
                }
                onevrest <- merge(onevrest, emobdt_jj[, c("level_use_join", "response"), with = FALSE],
                  by.y = "level_use_join", by.x = "term_1"
                )

                onevall <- as.data.table(suppressWarnings(emmeans::contrast(emob[emobdt_jj[, idx]], "eff", wts = emobdt_jj[, wts], parens = NULL, adjust = "none", data = dat, by = conditionvar)))
                onevall[, term_1 := gsub(" effect", "", contrast)]
                onevall[, contrast := gsub("effect", "vs. avg.all", contrast)]
                emobdt_jj[, level_use_join := level]
                if (length(confounded_covariates) > 0) {
                  for (ii in seq_along(confounded_covariates)) {
                    emobdt_jj[["level_use_join"]] <- paste0(emobdt_jj[["level_use_join"]], " ", emobdt_jj[[confounded_covariates[ii]]])
                    add_parenth <- onevall[grepl("\\(", term_1), term_1]
                    emobdt_jj[
                      level_use_join %in% gsub("\\(|\\)", "", add_parenth),
                      level_use_join := paste0("(", level_use_join, ")")
                    ]
                  }
                }
                onevall <- merge(onevall, emobdt_jj[, c("level_use_join", "response"), with = FALSE],
                  by.y = "level_use_join", by.x = "term_1"
                )
                emobdt_jj[, level_use_join := NULL]
                if (linkused == "identity") {
                  onevrest[, fold_change := response / (response - estimate)]
                  onevall[, fold_change := response / (response - estimate)]
                } else if (linkused == "log") {
                  onevrest[, fold_change := ratio]
                  onevall[, fold_change := ratio]
                  pw[, fold_change := ratio]
                } else {
                  onevrest[, fold_change := NA_real_]
                  onevall[, fold_change := NA_real_]
                  pw[, fold_change := NA_real_]
                }
                onevrest <- merge(onevrest,
                  stats_rest[level != "avg.rest"][, c("term_1", "counts", "propnz", "ncells", conditionvar), with = FALSE],
                  by = c("term_1", conditionvar),
                  sort = FALSE,
                  all.x = TRUE
                )

                onevrest <- merge(onevrest,
                  stats_rest[level == "avg.rest"][, c("term_1", "counts", "propnz", "ncells", conditionvar), with = FALSE],
                  by = c("term_1", conditionvar),
                  sort = FALSE,
                  all.x = TRUE,
                  suffixes = c("_1", "_2")
                )

                onevall <- merge(onevall,
                  stats_rest[level != "avg.rest"][, c("term_1", "counts", "propnz", "ncells", conditionvar), with = FALSE],
                  by = c("term_1", conditionvar),
                  sort = FALSE,
                  all.x = TRUE
                )

                onevall <- merge(onevall,
                  stats_all[, c("term_1", "counts", "propnz", "ncells", conditionvar), with = FALSE],
                  by = c("term_1", conditionvar),
                  suffixes = c("_1", "_2"),
                  all.x = TRUE,
                  sort = FALSE
                )

                onevall[, modelest_cpc_1 := response]
                onevrest[, modelest_cpc_1 := response]
                onevall[, modelest_cpc_2 := modelest_cpc_1 / fold_change]
                onevrest[, modelest_cpc_2 := modelest_cpc_1 / fold_change]
                rm_names <- c("response", "term_1")
                onevall[, (rm_names) := NULL]
                onevrest[, (rm_names) := NULL]

                onevrest_list[[term]] <- rbindlist(list(
                  onevrest_list[[term]],
                  onevrest[, `:=`(term = term)]
                ), use.names = TRUE, fill = TRUE)

                onevall_list[[term]] <- rbindlist(list(
                  onevall_list[[term]],
                  onevall[, `:=`(term = term)]
                ), use.names = TRUE, fill = TRUE)
              }
              ## Manually calculate fold change ratios for models with identity links (i.e., gaussian)
              if (linkused == "identity") {
                fc_pw <-
                  merge(emobdt_jj[, c("level", "response", "ky", "category"), with = FALSE],
                    emobdt_jj[, c("level", "response", "ky", "category"), with = FALSE],
                    by = "ky", allow.cartesian = TRUE, suffixes = c("_1", "_2")
                  )
                fc_pw <- fc_pw[category_1 < category_2]
                fc_pw[, fold_change := response_1 / response_2]
                fc_pw[, contrast := paste0(level_1, " / ", level_2)]
                pw[, level_1 := tstrsplit(contrast, " - ")[[1]]]
                pw[, level_2 := tstrsplit(contrast, " - ")[[2]]]
                pw <- merge(pw, fc_pw[, .(level_1, level_2, fold_change)],
                  by = c("level_1", "level_2")
                )

                rm_names <- c("level_1", "level_2")
                pw[, (rm_names) := NULL]
              } else if (linkused == "log") {
                pw[, fold_change := ratio]
              } else {
                pw[, fold_change := NA_real_]
              }
              ### Add summary stats to output
              pw[, term_1 := tstrsplit(contrast, "\\ [\\/|-]\\ ")[[1]]]
              pw[, term_2 := tstrsplit(contrast, "\\ [\\/|-]\\ ")[[2]]]
              pw <- merge(pw,
                stats_rest[level != "avg.rest"][, c("term_1", "counts", "propnz", "ncells", conditionvar), with = FALSE],
                by = c("term_1", conditionvar),
                sort = FALSE, all.x = TRUE
              )
              pwstats <- copy(stats_rest)
              setnames(pwstats, old = c("term_1"), new = c("term_2"))
              pw <- merge(pw,
                pwstats[level != "avg.rest"][, c("term_2", "counts", "propnz", "ncells", conditionvar), with = FALSE],
                by = c("term_2", conditionvar),
                sort = FALSE, all.x = TRUE,
                suffixes = c("_1", "_2")
              )
              pw <- merge(pw, emobdt_jj[, .(term_1 = level, modelest_cpc_1 = response)],
                by = c("term_1"),
                all.x = TRUE, sort = FALSE
              )
              pw <- merge(pw, emobdt_jj[, .(term_2 = level, modelest_cpc_2 = response)],
                by = c("term_2"),
                all.x = TRUE, sort = FALSE
              )
              pw[, `:=`(term_1 = NULL, term_2 = NULL)]
              pw_list[[term]] <- rbindlist(list(
                pw_list[[term]],
                pw[, `:=`(term = term)]
              ), use.names = TRUE, fill = TRUE)
              emobdt_jj[, ky := NULL]
              emm_list[[term]] <- rbindlist(list(
                emm_list[[term]],
                emobdt_jj[, `:=`(term = term)]
              ), use.names = TRUE, fill = TRUE)
            }
          }
        }
      }
    }
  }
  if (is.null(error_msg)) {
    if (!fitmethod %in% c("nebula::nebula", "spaMM::fitme")) {
      msumm <- .coef_summary(mod, original_formula, fitmethod)
    } else {
      msumm <- switch(fitmethod,
        "nebula::nebula" = nebula_summary(mod)$msumm,
        "spaMM::fitme" = spamm_summary(mod)$msumm
      )
    }
  } else {
    msumm <- data.table(msg = error_msg)
  }
  if (!missing(modtime)) {
    msumm <-
      tryCatch(
        {
          suppressMessages(
            cbind(msumm, rbind(modtime))[, `:=`(user.child = NULL, sys.child = NULL)]
          )
        },
        error = function(e) {
          return(msumm)
        }
      )
  }

  return(list(
    model_summary = msumm,
    pairwise_list = pw_list,
    onevrest_list = onevrest_list,
    onevall_list = onevall_list,
    emm_list = emm_list,
    terms = allterms
  ))
}
