library(testthat)
library(InSituDiff)
rm(list = ls())
data(minisle)
annot <- minisle$annot
rownames(annot) <- minisle$annot$cell_ID
counts <- minisle$counts
xy <- minisle$xy
norm <- sweep(counts, 1, Matrix::rowSums(counts)+10, "/")

# for running functions line by line:
if (FALSE) {
  mat = norm 
  xy = xy
  neighbors = NULL
  tissue = annot$tissue
  controlnames = c("control1", "control2") 
  k = 50
  radius = NULL
  ncontolneighborhoods = 1
  residtype = "diff"
  verbose = TRUE
}


#### test initialize function ------------------------------
obj <- suppressWarnings(initializeISD(
  mat = norm, # matrix of normalized expression
  xy = xy, # matrix of xy coordinates
  tissue = annot$tissue, # vector of cells' tissue IDs
  iscontrol = grepl("control", annot$tissue), 
  k = 50)) # define cellular neighborhoods as a cell's 50 nearest neighbors
test_that("initializeISD returns the expected results", {
  expect_true(all(dim(obj$neighbors) == nrow(norm)))
  expect_true(all(is.element(obj$controlmatches, rownames(norm)[grepl("control", annot$tissue)])))
})
#### test perturbation wrapper -----------------------------------


test_that("getPerturbation returns the expected results with linear diff", {
  selectedperturbations <- getPerturbations(
    x = norm, 
    obj = obj, 
    cells = rownames(counts)[annot$tissue == "SP19_1139"], # only look at one tissue
    genes = colnames(norm)[1:5],  
    residtype = "diff")
  expect_true(all(!is.na(selectedperturbations)))
  expect_true(all(dim(selectedperturbations) == c(sum(annot$tissue == "SP19_1139"), 5)))
})
test_that("getPerturbation returns the expected results with log2ratio diff", {
  selectedperturbations <- getPerturbations(
    x = norm, 
    obj = obj, 
    cells = rownames(counts)[annot$tissue == "SP19_1139"], # only look at one tissue
    genes = colnames(norm)[1:5],  
    residtype = "log2ratio")
  expect_true(all(!is.na(selectedperturbations)))
  expect_true(all(dim(selectedperturbations) == c(sum(annot$tissue == "SP19_1139"), 5)))
})


#### test subcomponents of perturbation calcs -------------------------

# test neighbor identification
test_that("getNeighbors returns the expected results", {
  set.seed(0)
  randominds <- sample(1:nrow(xy))
  tempneighbors <- getNeighbors(xy = xy[randominds, ], 
                                neighbors = NULL, 
                                tissue = annot$tissue, 
                                k = 50, 
                                radius = NULL, 
                                verbose = FALSE)
  # right number of neighbors:
  expect_true(all(Matrix::rowSums(tempneighbors != 0) == 50))
  # no neighbors from other tissues:
  expect_true(all(tempneighbors[annot$tissue == "control1", annot$tissue == "control2"] == 0) )
  expect_true(all(tempneighbors[annot$tissue == "control1", annot$tissue == "SP19_1139"] == 0) )
  expect_true(all(tempneighbors[annot$tissue == "control2", annot$tissue == "SP19_1139"] == 0) )
})

# test matchToControls:
test_that("matchToControls returns the expected results", {
  neighbors <- getNeighbors(xy = xy, 
                            neighbors = NULL, 
                            tissue = annot$tissue, 
                            k = 50, 
                            radius = NULL, 
                            verbose = FALSE) 
  neighbormat <- getNeighborhoodExpression(x = norm, 
                                           neighbors = neighbors,
                                           makedense = TRUE)
  nmatcolmeans <- Matrix::colMeans(neighbormat)
  neighbormat <- sweep(neighbormat, 2, pmax(nmatcolmeans, quantile(nmatcolmeans, 0.5)), "/") 
  controlmatches <- matchToControls(x = neighbormat, 
                                    tissue = annot$tissue, 
                                    iscontrol = is.element(annot$tissue, c("control1", "control2"))) 
  
  # matches are all from the right subset of controls:
  expect_true(all(annot$tissue[controlmatches] != "SP19_1139"))
  expect_true(all(annot$tissue[controlmatches[annot$tissue == "control1"]] == "control2"))
  expect_true(all(annot$tissue[controlmatches[annot$tissue == "control2"]] == "control1"))
})



# test module search:
modules <- buildGeneModules(
  x = norm, obj = obj, 
  genes = NULL,  # subset of genes to use; enter NULL to use all, or "highlyperturbed" to take just the most perturbed genes ("highlyperturbed" is our default and our recommendation for full studies; here we use NULL because our mini dataset has only 30 genes.)
  resolution = 0.02, # control leiden clustering
  corthresh = 0.1, # control the adjacency network used by leiden clustering
  min_module_cor = 0.1, # throw out modules with average cor below this
  subsetsize = 1e5, 
  eps = 0.01, 
  residtype = "log2ratio") 
test_that("buildGeneModules returns the expected results", {
  expect_true(is.list(modules))
})

# test scoring of module perturbations:
test_that("getPerturbation returns the expected results if given modules", {
  selectedperturbations <- getPerturbations(
    x = norm, 
    obj = obj, 
    cells = rownames(counts)[annot$tissue == "SP19_1139"], # only look at one tissue
    genes = modules[1:2],  
    residtype = "log2ratio")
  expect_true(all(!is.na(selectedperturbations)))
  expect_true(all(dim(selectedperturbations) == c(sum(annot$tissue == "SP19_1139"), 2)))
})

# test spatial clustering:
test_that("buildGeneModules returns the expected results", {
  res <- clusterPerturbations(x = norm, 
                              obj = obj, 
                              cells = NULL, 
                              genes = colnames(norm),
                              residtype = "log2ratio", 
                              nclust = 6,
                              eps = 0.1)
  expect_true(is.vector(res$clust))
  expect_true(is.matrix(res$means))
})