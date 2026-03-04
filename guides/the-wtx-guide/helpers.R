# --- Helper: splitScale ---
#' @title Find Medoids for Cluster Labeling
#' @description Finds the medoids (most central points via PAM clustering) of groups 
#'   defined in coords_df, allowing k > 1 per group for sub-cluster labeling.
#' @param coords_df data.frame containing coordinates (x_coords, y_coords) and the label column.
#' @param k_df data.frame with columns 'group', 'group_label', and 'k'.
#' @param upper_limit numeric. Max number of cells to sample per group for PAM clustering (max 65536).
#' @return A data.frame ready for plotting labels, containing the 'index', 'x_coords', 'y_coords', 'group', and 'group_label' of all medoid points.
findMedoids <- function(coords_df, k_df, upper_limit=1e4){
  require(cluster)
  require(purrr)
  
  label_by <- setdiff(colnames(coords_df), c("Run_Tissue_name", "x_coords", "y_coords"))
  if (length(label_by) != 1) stop("Could not unambiguously determine the label column.")
  
  if(upper_limit > 65536) upper_limit <- 65536
  if(!inherits(k_df, 'data.frame')) stop("k_df must be a data.frame.")
  
  if(nrow(k_df)==0){
    groups <- levels(as.factor(coords_df[[label_by]]))
    k_df <- data.frame(
      group = groups,
      group_label = groups,
      k = c(rep(1, length(groups)))
    )
  } else if(!all(c("group", "group_label", "k") %in% colnames(k_df))){
    stop(paste0("k_df needs columns: group, group_label, k."))
  }
  
  coords_df$.index <- 1:nrow(coords_df)
  
  medoids_list <- 1:nrow(k_df) %>%
    purrr::map(function(i){
      
      group_name <- k_df[i,'group']
      group_label <- k_df[i,'group_label']
      k_val <- k_df[i,'k']
      
      coords_subset <- coords_df[coords_df[,label_by] == group_name,]
      if(nrow(coords_subset) == 0) return(NULL)
      
      if(nrow(coords_subset) == 1){
        # Edge case: single cell
        return(data.frame(index = coords_subset[1,'.index'],
                          group = group_name, x_coords = coords_subset[1,'x_coords'],
                          y_coords = coords_subset[1,'y_coords'], group_label = group_label))
      }
      
      # Reduce size if needed
      if(nrow(coords_subset) > upper_limit){
        set.seed(98103)
        warning(paste0("Subset needed for group ", group_name, ". Sampling ", upper_limit, " cells."), call. = FALSE)
        coords_subset <- coords_subset[sample(x=1:nrow(coords_subset), size = upper_limit, replace=FALSE), ] 
      }
      
      # Lower-level function call (PAM)
      meds <- pam(coords_subset[,c('x_coords', 'y_coords')], k=k_val)
      medoids_cells <- coords_subset[meds$id.med, ]
      
      return(data.frame(index = medoids_cells$.index,
                        group = rep(group_name, nrow(medoids_cells)),
                        x_coords = medoids_cells$x_coords,
                        y_coords = medoids_cells$y_coords,
                        group_label = rep(group_label, nrow(medoids_cells)))
      )
    })
  
  medoids_df <- do.call(rbind, medoids_list)
  
  return(medoids_df)
}

splitScale <- function(df, n, scale_colors){
  if(n==1){
    df$color <- scale_colors[1]
  } else {
    for(i in seq(from=2, to=n)){
      new_x <- df[nrow(df),'xstart'] + (df[nrow(df),'xend'] - df[nrow(df),'xstart']) / 2
      df <- rbind(df, df[nrow(df),])
      df[nrow(df)-1, 'xend'] <- new_x
      df[nrow(df), 'xstart'] <- new_x
    }
    df$color=rep(scale_colors, nrow(df))[1:nrow(df)]
  }
  
  return(df)
}

addSimpleScale <- function(g, location=c(9, 6.5), padding_proportion=0.1, width=1, height=0.05, n=1,
                           scale_colors=c("black"), display_units = "mm",
                           label_nudge_x=0, label_nudge_y= -0.1, label_size=3){
  
  ### ###########
  ### Check input
  ### ###########
  
  if(!inherits(g, "ggplot")){ stop("g must be a ggplot object.") }
  if(is.numeric(location)){ if(length(location) != 2) stop('Numeric locations need to be of the form c(x, y) only.') } 
  else if(is.character(location)){ 
    if(length(location) != 1) stop("Expecting a single chacter for location.")
    acceptable_locations <- c("bl")
    if(!location %in% acceptable_locations) stop(paste0("Only the following acceptable locations are allowed currently: ", paste(acceptable_locations, collapse = ", "), "."))
    if(!inherits(padding_proportion, "numeric") | padding_proportion <0 | padding_proportion > 1) stop("proportions are between 0 and 1 (inclusive).")
  } else { stop("Was expecting either a dublet showing the x and y location of scale bar or an acceptable abbreviation of the location.") }
  if(!inherits(width, "numeric") | !inherits(height, "numeric")) stop('width and height need to be numeric')
  if(!is.numeric(n) | n %% 1 != 0) n <- as.integer(n)
  if(n>1 & length(scale_colors)!=2) stop("Need 2 colors if n > 1")
  if(!display_units %in% c("mm", "µm")) stop("display units supported are mm and µm only.")
  
  ### #######
  ### Process
  ### #######

  stopifnot('geom_point' %in% names(g$layers))
  plot_data <- g$layers$geom_point$data

  if (inherits(plot_data, "waiver")) plot_data <- g$data
  if (is.null(plot_data) || inherits(plot_data, "waiver")) {
    warning("Could not find data in ggplot object; skipping scale bar addition.", call. = FALSE)
    return(g) 
  }

  x <- plot_data[['x_coords']] 
  y <- plot_data[['y_coords']]
  

  if(is.numeric(location)){
    df <- data.frame(xstart=location[1], xend=location[1] + width, ystart=location[2], yend=location[2] + height)
    text_x <- df$xstart[1] + (df$xend[1] - df$xstart[1])/2 + label_nudge_x
    text_y <- df$ystart[1] + label_nudge_y
  } else if(location=="bl"){
    df <- data.frame(xstart=min(x)*(1+padding_proportion), xend=min(x)*(1+padding_proportion) + width, 
                     ystart=min(y)*(1+padding_proportion), yend=min(y)*(1+padding_proportion) + height)
    text_x <- df$xstart[1] + (df$xend[1] - df$xstart[1])/2 + label_nudge_x
    text_y <- df$ystart[1]*(1-padding_proportion) + label_nudge_y
  }
  
  df <- splitScale(df, n, scale_colors)
  fill_color_vec <- scale_colors
  names(fill_color_vec) <- scale_colors
  text_label <- ifelse(display_units=="mm", paste0(width, " mm"), paste0(width*1e3, " µm"))
  
  p <- g + geom_rect(data=df, 
                     aes(xmin=xstart, xmax=xend, ymin=ystart, ymax=yend, fill=color),
                     inherit.aes = FALSE) + 
    scale_fill_manual(values=c(fill_color_vec), guide="none")
  
  p <- p + annotate(geom = "text",
                    x = text_x, 
                    y = text_y,
                    label = text_label, 
                    size = label_size,
                    fontface = 'plain')
  
  return(p)
}

setGeneric(name="plotDots", signature = c("input"),
           function(input, ...){
             standardGeneric("plotDots")
           })

setMethod(f="plotDots",
          signature="anndata._core.anndata.AnnData",
          definition = function(input, 
                                obsm_key = "spatial",
                                plot_global = TRUE,
                                color_by = NULL,
                                color_is_numeric = NULL, 
                                color_layer = NULL,
                                facet_by_group = FALSE,
                                additional_ggplot_layers = list(),
                                additional_plot_parameters = list()
          ){

            require(reticulate)
            require(utils)

            obs_keys <- reticulate::py_to_r(input$obs_keys())
            obsm_keys <- reticulate::py_to_r(input$obsm_keys())
            stopifnot(obsm_key %in% obsm_keys)
            is_spatial <-  obsm_key == 'spatial'
            
            base_obs_cols <- c("Run_Tissue_name")
            is_in_obs <- !is.null(color_by) && color_by %in% obs_keys
            if (is_in_obs) {
              obs_cols_to_get <- unique(c(base_obs_cols, color_by))
            } else {
              obs_cols_to_get <- base_obs_cols
            }
            
            missing_obs <- setdiff(obs_cols_to_get, obs_keys)
            if (length(missing_obs) > 0) stop(paste("Missing columns in adata.obs:", paste(missing_obs, collapse=", ")))
            
            df <- reticulate::py_to_r(input$obs[obs_cols_to_get])
            if(!"Run_Tissue_name" %in% colnames(df)) {
              df$Run_Tissue_name <- "S0"
            }
            
            coords_np <- input$obsm[obsm_key]
            coord_data_matrix <- reticulate::py_to_r(coords_np)
            if (ncol(coord_data_matrix) < 2) stop(paste0("obsm_key '", 
                                  obsm_key, "' has fewer than 2 dimensions."))
            coord_data <- as.data.frame(coord_data_matrix[, 1:2])
            colnames(coord_data) <- c("x_coords", "y_coords")
            
            if (nrow(df) != nrow(coord_data)) stop("Row count mismatch between adata.obs and adata.obsm.")
            df <- cbind(df, coord_data)
            
            if (!is_in_obs && !is.null(color_by)) {
              if (!color_by %in% reticulate::py_to_r(input$var_names$to_list())) stop(paste0("Gene '", color_by, "' not found in adata.var_names."))
              all_exprs_np <- input$obs_vector(color_by, layer = color_layer)
              all_exprs_r <- reticulate::py_to_r(all_exprs_np)
              if(length(all_exprs_r) != nrow(df)) stop("Row count mismatch between df and expression data.")
              df[[color_by]] <- all_exprs_r
            }
          
            plotDots(input = df,
                     coords_key = obsm_key,
                     is_spatial = is_spatial,
                     plot_global = plot_global,
                     color_by = color_by,
                     color_is_numeric = color_is_numeric,
                     color_layer = color_layer,
                     facet_by_group = facet_by_group,
                     additional_ggplot_layers = additional_ggplot_layers,
                     additional_plot_parameters = additional_plot_parameters)
          })

setMethod(f="plotDots",
          signature="data.frame",
          definition = function(input, 
                                coords_key = NULL,
                                is_spatial = NULL,
                                plot_global = TRUE,
                                color_by = NULL,
                                color_is_numeric = NULL, 
                                color_layer = NULL, 
                                facet_by_group = FALSE,
                                additional_ggplot_layers = NULL,
                                additional_plot_parameters = list()
          ){
            require(ggplot2)
            require(purrr)
            require(utils) 
            require(ggrepel)
            
            df <- input
            
            if(!all(c('x_coords', 'y_coords') %in% colnames(df))){
              stop("input needs to have columns x_coords and y_coords at a minimum.")
            }
            
            default_params <- list(directory = getwd(), prefix = "plotDots_", width = 4, height = 4, 
                                   fileType = "png", dpi = 100,
                                   geom_point_params = list(size = 0.5, alpha = 1),
                                   geom_label_params = list(size = 5, alpha = 1, min.segment.length = 0.1,
                                                            fontface = 'bold', max.overlaps = 40,
                                                            box.padding = unit(0.35, "lines"),
                                                            point.padding = unit(0.5, "lines"),
                                                            segment.color = 'grey50'))
            params <- utils::modifyList(default_params, additional_plot_parameters)
            
            # Check for Multi-Slide Execution
            unique_slides <- unique(df$Run_Tissue_name)
            
            if (is_spatial && length(unique_slides) > 1) {

              message(paste("Initiating multi-slide processing for", length(unique_slides), "slides..."))
              
              plot_list <- unique_slides %>% 
                purrr::map(~ {
                  df_subset <- df[df$Run_Tissue_name == .x, ]
                  
                  plotDots(input = df_subset,
                           coords_key = coords_key,
                           is_spatial = is_spatial, color_by = color_by,
                           color_is_numeric = color_is_numeric, color_layer = color_layer,
                           facet_by_group = facet_by_group, 
                           additional_ggplot_layers = additional_ggplot_layers,
                           additional_plot_parameters = additional_plot_parameters)
                })
              return(invisible(plot_list))
              
            } else {
              
              if (is.null(color_by)) {
                color_col <- NULL
                color_is_numeric <- FALSE
              } else {
                color_col <- color_by
                if (is.null(color_is_numeric)) color_is_numeric <- is.numeric(df[[color_col]])
                if (color_is_numeric) df <- df[order(df[[color_col]], decreasing = FALSE), ]
              }
              slide_id <- unique(df$Run_Tissue_name)
              

              if('labels_on_plot' %in% names(params) && !color_is_numeric){
                medoids_df <- findMedoids(df, k_df=params$labels_on_plot)
                use_medoids <- TRUE
              } else {
                use_medoids <- FALSE
              }
              

              if(plot_global){
                if(is.null(color_col)){
                  p_base <- ggplot(data = df, 
                                   aes(x = x_coords, y = y_coords))                  
                } else {
                  p_base <- ggplot(data = df, 
                                   aes(x = x_coords, y = y_coords, 
                                       color = .data[[color_col]]))
                }


                p_global <- p_base + rlang::exec(ggplot2::geom_point, !!!params$geom_point_params) 

                if (!is.null(color_col) && color_is_numeric) {
                  p_global <- p_global + scale_color_viridis_c(limits = range(df[[color_col]]))
                }
                
                p_global <- p_global + additional_ggplot_layers
                

                if (is_spatial && 'scale_bar_params' %in% names(params)) {
                  scale_args <- c(list(g = p_global), params$scale_bar_params)
                  p_global <- tryCatch({
                    rlang::exec(addSimpleScale, !!!scale_args)
                  }, error = function(e) {
                    warning(paste("Scale bar failed with error:", e$message), call. = FALSE)
                    p_global
                  })
                }

                if(use_medoids){

                    fixed_aes <- ggplot2::aes(x = x_coords, y = y_coords, 
                                              label = as.character(group_label),
                                              color = as.character(group))
                    
                    label_args <- c(
                      list(data = medoids_df,
                           mapping = fixed_aes), 
                      params$geom_label_params 
                    )
                    
                    p_global <- p_global + rlang::exec(ggrepel::geom_label_repel, !!!label_args)
                }
                  

                if(!is.null(params$fileType)){
                  filename <- paste0(params$prefix, "__", coords_key, "__", slide_id, "__plot.", params$fileType)
                  ggsave(filename = file.path(params$directory, filename), plot = p_global, 
                         width = params$width, height = params$height, device = params$fileType, dpi = params$dpi)                  
                }

              }
              

              if (facet_by_group && !is.null(color_col) && !color_is_numeric) {
                
                p_facet_list <- levels(as.factor(df[[color_col]])) %>%
                  purrr::map(~ {
                    group_name <- .x
                    df_focus <- df[df[[color_col]] == group_name, ]
                    
                    p_facet <- ggplot() + 

                      geom_point(data = df, aes(x = x_coords, y = y_coords),
                                 color = 'grey80', size = 0.5, alpha = 0.1) +
                      

                      rlang::exec(ggplot2::geom_point, data = df_focus, 
                                  aes(x = x_coords, y = y_coords, color = .data[[color_col]]), 
                                  !!!params$geom_point_params) +
                      
                      ggplot2::scale_color_hue() + 
                      
                      additional_ggplot_layers + 
                      
                      theme(legend.position = "none") + labs(title = group_name)
                    
                    if (is_spatial && 'scale_bar_params' %in% names(params)) {
                      p_facet <- rlang::exec(addSimpleScale, g = p_facet, !!!params$scale_bar_params)
                    }
                    
                    if(!is.null(params$fileType)){
                      facet_filename <- paste0(params$prefix, "__", coords_key, "__", slide_id, "__facet_", group_name, ".", params$fileType)
                      ggsave(filename = file.path(params$directory, facet_filename), plot = p_facet, 
                             width = params$width, height = params$height, device = params$fileType, dpi=params$dpi)
                    }
                    return(p_facet)
                  })
                

                return(invisible(c(if(exists("p_global")) list(Global = p_global) else NULL, 
                                   p_facet_list)))
              }
              
              return(invisible(if(exists("p_global")) p_global else NULL))
            }
          })

