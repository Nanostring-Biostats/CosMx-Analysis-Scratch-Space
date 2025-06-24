test_that("multiplication works", {
  expect_equal(2 * 2, 4)
})

test_that("pre_de outputs", {
  require(smiDE)
  require(data.table)
  datadir<-system.file("extdata", package="smiDE")
  sem <- readRDS(paste0(datadir, "/small_nsclc.rds"))
  totalcounts <- sem@meta.data[["totalcounts"]]
  sem <- Seurat::SetAssayData(sem
                              ,"data"
                              ,sem[["RNA"]]@counts %*% Matrix::Diagonal(x=mean(totalcounts)/totalcounts, names=colnames(sem))
  )
  metainfo <- data.table(sem@meta.data)                           
 
  ref_celltypes <- c("macrophage", "fibroblast") 
    pre_de_obj <- 
      pre_de(counts = sem[["RNA"]]@counts
             ,normalized_data = sem[["RNA"]]@data
             ,metadata = metainfo
             ,cell_type_metadata_colname = "cell_type"
             ,split_neighbors_by_colname = "tissue"
             ,mm_radius = 0.05
             ,ref_celltype = ref_celltypes
             ,sdimx_colname = "sdimx"
             ,sdimy_colname = "sdimy"
             ,weight_colname = "weight"
             ,aggregation = "sum"
             ,verbose=TRUE
             ,adjacencies_only = FALSE
      )
      
 
    
    all_cell_ids <-  metainfo[["cell_ID"]] 
    cell_ids <- metainfo[cell_type %in% ref_celltypes][["cell_ID"]] 
    expect_true(all(c("neighbor_expr_byct", "adjacency_mat", "ct_matrix", "adjacency_counts_by_ct", "ref_celltype") 
                    %in% names(pre_de_obj$nblist)))

    expect_true(all(all_cell_ids %in% rownames(pre_de_obj$nblist$adjacency_mat)))
    expect_true(all(all_cell_ids %in% colnames(pre_de_obj$nblist$adjacency_mat)))
    expect_true(length(all_cell_ids) == dim(pre_de_obj$nblist$adjacency_mat)[1])
    expect_true(dim(pre_de_obj$nblist$adjacency_mat)[1] == dim(pre_de_obj$nblist$adjacency_mat)[2])
  
    ### check adjacency_counts_by_ct 
    expect_true(all(metainfo[["cell_type"]] %in% colnames(pre_de_obj$nblist$adjacency_counts_by_ct))) 
    expect_true(all(all_cell_ids %in% pre_de_obj$nblist$adjacency_counts_by_ct[["cell_ID"]]))
   
    ### check ct_matrix 
    expect_true(all(metainfo[["cell_type"]] %in% colnames(pre_de_obj$nblist$ct_matrix))) 
    expect_true(all(all_cell_ids %in% rownames(pre_de_obj$nblist$ct_matrix)))
    expect_true(all(Matrix::rowSums(pre_de_obj$nblist$ct_matrix)==1))
    
    ### check that neighborhood counts objects are there for each ct + "allct" and "otherct", 
    ### and have the correct dimensions and dimnames 
    expect_true(all(cell_ids %in% colnames(pre_de_obj$nblist$neighbor_expr_byct[["otherct"]])))
    expect_true(all(rownames(sem[["RNA"]]@counts) %in% rownames(pre_de_obj$nblist$neighbor_expr_byct[["otherct"]])))
    expect_true(all(rownames(sem[["RNA"]]@counts) %in% rownames(pre_de_obj$nblist$neighbor_expr_byct[["allct"]])))
  
    uniq_ct <- unique(metainfo[["cell_type"]])
    
    for(ct in uniq_ct){
      ctname <- gsub("\\ |\\-", "_", ct)
      expect_true(ctname %in% 
                    names(pre_de_obj[["nblist"]][["neighbor_expr_byct"]]))
      expect_equal(dimnames(pre_de_obj[["nblist"]][["neighbor_expr_byct"]][[ctname]])
                  ,dimnames(pre_de_obj[["nblist"]][["neighbor_expr_byct"]][["otherct"]])
                  )
    }
    
     
})





