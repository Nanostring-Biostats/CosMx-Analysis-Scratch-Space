test_that("smi_de and results", {
  require(smiDE)
  require(data.table)
  datadir <- system.file("extdata", package="smiDE")
  sem <- readRDS(paste0(datadir, "/small_nsclc.rds"))
  totalcounts <- sem@meta.data[["totalcounts"]]
  sem <- Seurat::SetAssayData(sem
                             ,"data"
                             ,sem[["RNA"]]@counts %*% Matrix::Diagonal(x=mean(totalcounts)/totalcounts, names=colnames(sem))
                             )
  metainfo <- data.table(sem@meta.data)                           
  
  pre_de_obj <- 
  pre_de(counts = sem[["RNA"]]@counts
         ,normalized_data = sem[["RNA"]]@data
         ,metadata = metainfo
         ,cell_type_metadata_colname = "cell_type"
         ,split_neighbors_by_colname = "tissue"
         ,mm_radius = 0.05
         ,ref_celltype = c("macrophage", "fibroblast")
         ,sdimx_colname = "sdimx"
         ,sdimy_colname = "sdimy"
         ,weight_colname = "weight"
         ,aggregation = "sum"
         ,verbose=TRUE
         ,adjacencies_only = FALSE
  )
  
  pre_de_obj_barebones <- 
  pre_de(metadata = metainfo
         ,adjacencies_only = TRUE
         ,cell_type_metadata_colname = "cell_type"
         ,split_neighbors_by_colname = "tissue"
         ,mm_radius = 0.05
         ,sdimx_colname = "sdimx"
         ,sdimy_colname = "sdimy"
         ,verbose=TRUE
  )
 
  fibroblast_and_macrophage_cells <- metainfo[cell_type %in% c(pre_de_obj$nblist$ref_celltype),cell_ID]
  
  ## neg binomial , FE model
  expect_warning(
  de_nb_fe <-       
    smi_de(assay_matrix = sem[["RNA"]]@counts[,fibroblast_and_macrophage_cells]
           ,metadata = metainfo[cell_ID %in% fibroblast_and_macrophage_cells]
           ,formula = ~RankNorm(otherct_expr) + niche + tissue + offset(log(totalcounts)) 
           ,neighborhood_counts = pre_de_obj$nblist
           ,family="nbinom2"
           ,targets=rownames(sem)[1:3]
    ) 
    
  )
  
  scalefactors_for_neighbor_expression <- 
    metainfo[match(colnames(sem), cell_ID),mean(totalcounts)/totalcounts] 
  expect_error( ### didn't name the scalefactors with the cell ids 
    de_nb_fe_prede_onthefly <-       
      smi_de(assay_matrix = sem[["RNA"]]@counts
             ,metadata = metainfo[cell_ID %in% fibroblast_and_macrophage_cells]
             ,formula = ~RankNorm(otherct_expr) + niche + tissue + offset(log(totalcounts)) 
             ,pre_de_obj = pre_de_obj_barebones
             ,neighbor_expr_cell_type_metadata_colname = "cell_type"
             ,neighbor_expr_overlap_weight_colname = NULL
             ,neighbor_expr_overlap_agg ="sum"
             ,neighbor_expr_totalcount_normalize = TRUE
             ,neighbor_expr_totalcount_scalefactor = scalefactors_for_neighbor_expression
             ,family="nbinom2"
             ,targets=rownames(sem)[1:3]
      ) 
      
  ) 
  names(scalefactors_for_neighbor_expression) <- colnames(sem)
  de_nb_fe_prede_onthefly <-       
    smi_de(assay_matrix = sem[["RNA"]]@counts
           ,metadata = metainfo[cell_ID %in% fibroblast_and_macrophage_cells]
           ,formula = ~RankNorm(otherct_expr) + niche + tissue + offset(log(totalcounts)) 
           ,pre_de_obj = pre_de_obj_barebones
           ,neighbor_expr_cell_type_metadata_colname = "cell_type"
           ,neighbor_expr_overlap_weight_colname = NULL
           ,neighbor_expr_overlap_agg ="sum"
           ,neighbor_expr_totalcount_normalize = TRUE
           ,neighbor_expr_totalcount_scalefactor = scalefactors_for_neighbor_expression
           ,family="nbinom2"
           ,targets=rownames(sem)[1:3]
    ) 
  
  #### compare pre-calcualted pre_de_object vs. 'on-the-fly'
  exclude_cols <- c("user.self", "sys.self", "elapsed")
  expect_equal(results(de_nb_fe, comparisons = "model_summary")[[1]][,-c(exclude_cols),with=FALSE]
               ,results(de_nb_fe_prede_onthefly, comparisons = "model_summary")[[1]][,-c(exclude_cols),with=FALSE]
  )
  for(comp in c("pairwise", "one.vs.rest", "one.vs.all")){
    expect_equal(results(de_nb_fe, comparisons = comp)
                 ,results(de_nb_fe_prede_onthefly, comparisons = comp)
                 )
      
  } 
  ## neg binomial , RE model
  expect_warning(
  de_nb_re <-       
    smi_de(assay_matrix = sem[["RNA"]]@counts[,fibroblast_and_macrophage_cells]
           ,metadata = metainfo[cell_ID %in% fibroblast_and_macrophage_cells]
           ,formula = ~RankNorm(otherct_expr) + niche + (1|tissue) + offset(log(totalcounts)) 
           ,neighborhood_counts = pre_de_obj$nblist
           ,family="nbinom2"
           ,targets=rownames(sem)[1:3]
    ) 
  ) 
  
  ## gaussian , FE model
  expect_warning(
  de_gaussian_fe <-       
    smi_de(assay_matrix = sem[["RNA"]]@data[,fibroblast_and_macrophage_cells]
           ,metadata = metainfo[cell_ID %in% fibroblast_and_macrophage_cells]
           ,formula = ~RankNorm(otherct_expr) + niche + tissue 
           ,neighborhood_counts = pre_de_obj$nblist
           ,family="gaussian"
           ,targets=rownames(sem)[1:3]
    ) 
  )
  
  ## gaussian , RE model
  expect_warning(
  de_gaussian_re <-       
    smi_de(assay_matrix = sem[["RNA"]]@data[,fibroblast_and_macrophage_cells]
           ,metadata = metainfo[cell_ID %in% fibroblast_and_macrophage_cells]
           ,formula = ~RankNorm(otherct_expr) + niche + (1|tissue) 
           ,neighborhood_counts = pre_de_obj$nblist
           ,family="gaussian"
           ,targets=rownames(sem)[1:3]
    ) 
  ) 
  ## poisson , FE model
  expect_warning(
  de_poisson_fe <-       
    smi_de(assay_matrix = sem[["RNA"]]@counts[,fibroblast_and_macrophage_cells]
           ,metadata = metainfo[cell_ID %in% fibroblast_and_macrophage_cells]
           ,formula = ~RankNorm(otherct_expr) + niche + tissue 
           ,neighborhood_counts = pre_de_obj$nblist
           ,family="poisson"
           ,targets=rownames(sem)[1:3]
    ) 
  ) 
  
  ## poisson , RE model 
  expect_warning(
  de_poisson_re <-       
    smi_de(assay_matrix = sem[["RNA"]]@counts[,fibroblast_and_macrophage_cells]
           ,metadata = metainfo[cell_ID %in% fibroblast_and_macrophage_cells]
           ,formula = ~RankNorm(otherct_expr) + niche + (1|tissue) 
           ,neighborhood_counts = pre_de_obj$nblist
           ,family="poisson"
           ,targets=rownames(sem)[1:3]
    ) 
    
  )

  ## gaussian , spatial RE model  spaMM
  expect_warning(
  de_gaussian_sre_spamm <-       
    smi_de(assay_matrix = sem[["RNA"]]@counts[,fibroblast_and_macrophage_cells]
           ,metadata = metainfo[cell_ID %in% fibroblast_and_macrophage_cells]
           ,formula = ~RankNorm(otherct_expr) + niche + (1|tissue) 
           ,neighborhood_counts = pre_de_obj$nblist
           ,family="gaussian"
           ,targets=rownames(sem)[1:3]
           ,spatial_model = list(
             name = "GP_Matern"
             ,k_prop_n = 0.001
             ,x_coord_col = "sdimx"
             ,y_coord_col = "sdimy"
             ,split_neighbors_by_colname = "Run_Tissue_name"
             ,spatial_random_effect = ~Matern(1 | sdimx_cluster + sdimy_cluster %in% Run_Tissue_name)
           )
    ) 
  )
  
  ## nb , spatial RE model spaMM
  expect_warning(
  de_nb_sre_spamm <-       
    smi_de(assay_matrix = sem[["RNA"]]@counts[,fibroblast_and_macrophage_cells]
           ,metadata = metainfo[cell_ID %in% fibroblast_and_macrophage_cells]
           ,formula = ~RankNorm(otherct_expr) + niche + (1|tissue)  + offset(log(totalcounts))
           ,neighborhood_counts = pre_de_obj$nblist
           ,family="nbinom2"
           ,targets=rownames(sem)[1:3]
           ,spatial_model = list(
             name = "GP_Matern"
             ,k_prop_n = 0.001
             ,x_coord_col = "sdimx"
             ,y_coord_col = "sdimy"
             ,split_neighbors_by_colname = "Run_Tissue_name"
             ,spatial_random_effect = ~Matern(1 | sdimx_cluster + sdimy_cluster %in% Run_Tissue_name)
           )
    ) 
  )
 
  if(FALSE){
    de_nb_sre_inla <- 
      smiDE::smi_de(assay_matrix = sem[["RNA"]]@counts[,fibroblast_and_macrophage_cells]
                    ,metadata = metainfo[cell_ID %in% fibroblast_and_macrophage_cells]
                    ,formula = ~RankNorm(otherct_expr) + niche + tissue + offset(log(totalcounts))
                    ,neighborhood_counts = pre_de_obj$nblist
                    ,groupVar = "niche"
                    ,family="nbinom2"
                    ,targets = rownames(sem)[1:3]
                    ,spatial_model = list(name="GP_INLA", 
                                         quantiles = c(0.025/3, 0.5, 1-0.025/3)
                    )
                    ,nCores=1
      )
    
    de_gaussian_sre_inla <- 
      smiDE::smi_de(assay_matrix = sem[["RNA"]]@counts[,fibroblast_and_macrophage_cells]
                    ,metadata = metainfo[cell_ID %in% fibroblast_and_macrophage_cells]
                    ,formula = ~RankNorm(otherct_expr) + niche + tissue 
                    ,neighborhood_counts = pre_de_obj$nblist
                    ,groupVar = "niche"
                    ,family="nbinom2"
                    ,targets = rownames(sem)[1:3]
                    ,spatial_model = list(name="GP_INLA", 
                                         quantiles = c(0.025/3, 0.5, 1-0.025/3)
                    )
                    ,nCores=1
      )
      
  } 
    

  ### Check that counts, cells, and other results stats align as expected and 
  ##  output is there for all expected contrasts. 
  ## i.e., choose(9, 2) pairwise contrasts for a categorical term with 9 levels
  ## i.e., 9 contrasts for one.vs.rest / one.vs.all  categorical term with 9 levels.
  ncat_niche <- metainfo[,uniqueN(niche)] 
  numcells_niche <- metainfo[,.N,by=.(niche)]
  num_genes <- 3 
   
  for(comp in c("pairwise", "one.vs.rest", "one.vs.all")){
    if(comp=="pairwise"){
      num_contrasts <- choose(ncat_niche, 2)
    }
    if(comp %in% c("one.vs.rest", "one.vs.all")){
      num_contrasts <- ncat_niche
    } 
    compare_gaussian <-  
    merge(results(de_gaussian_fe, comp)[[1]]
          ,results(de_gaussian_re, comp)[[1]]
          ,by=c("contrast", "term", "target")
          ,suffixes=c("_fe", "_re"))
    
    compare_nb <-  
    merge(results(de_nb_fe, comp)[[1]]
          ,results(de_nb_re, comp)[[1]]
          ,by=c("contrast", "term", "target")
          ,suffixes=c("_fe", "_re"))
   
    compare_nb_gaussian <-  
    merge(results(de_nb_fe, comp)[[1]]
          ,results(de_gaussian_fe, comp)[[1]]
          ,by=c("contrast", "term", "target")
          ,suffixes=c("_nb", "_gaussian"))
    
    
    expect_true(compare_gaussian[term=="niche",uniqueN(contrast)]==num_contrasts)
    expect_true(nrow(compare_gaussian[term=="niche"]) == num_contrasts*num_genes)
    
    expect_true(compare_nb[term=="niche",uniqueN(contrast)]==num_contrasts)
    expect_true(nrow(compare_nb[term=="niche"]) == num_contrasts*num_genes)
    
  
    expect_true(compare_gaussian[,all.equal(counts_1_re, counts_1_fe)])
    expect_true(compare_gaussian[,all.equal(ncells_1_re, ncells_1_fe)])
    expect_true(compare_gaussian[,all.equal(counts_2_re, counts_2_fe)])
    expect_true(compare_gaussian[,all.equal(ncells_2_re, ncells_2_fe)])
    expect_true(compare_gaussian[,cor(estimate_fe, estimate_re)] > 0.95)
    
    expect_true(compare_nb[term=="niche",uniqueN(contrast)]==num_contrasts)
    expect_true(nrow(compare_nb[term=="niche"]) == num_contrasts*num_genes)
    
  
    expect_true(compare_nb[,all.equal(counts_1_re, counts_1_fe)])
    expect_true(compare_nb[,all.equal(ncells_1_re, ncells_1_fe)])
    expect_true(compare_nb[,all.equal(counts_2_re, counts_2_fe)])
    expect_true(compare_nb[,all.equal(ncells_2_re, ncells_2_fe)])
    expect_true(compare_nb[,cor(ratio_fe, ratio_re)] > 0.95)
  
    expect_true(compare_nb_gaussian[,all.equal(ncells_1_gaussian, ncells_1_nb)])
    expect_true(compare_nb_gaussian[,all.equal(ncells_2_gaussian, ncells_2_nb)])
    
      
  } 
  
})


