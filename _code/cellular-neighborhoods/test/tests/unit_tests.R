library(testthat)
rm(list = ls())

data(cosmx_kidney)
annot <- cosmx_kidney$annot
rownames(annot) <- annot$cell_ID
counts <- cosmx_kidney$counts
celltype <- as.factor(cosmx_kidney$annot$celltype)
xy <- cosmx_kidney$xy
neighbors <- radiusBasedGraph(x = xy[, 1], y = xy[, 2], R = 0.05, subset=1)


#### test neighbor definition: -----------------------

neighbors <- radiusBasedGraph(x = xy[, 1], y = xy[, 2], R = 0.05, subset=1)
neighbors <- nearestNeighborGraph(x = xy[, 1], y = xy[, 2], N = 50, subset=1)
# (should be sparse matrices with dim = nrow(xy).)


#### test neighbor summing math: --------------------

test_that("neighbor math has the right logic", {
  
  expect_equal(neighbor_tabulate(annot$celltype, neighbors)[1:5, 1:3],
               neighbor_tabulate(annot$celltype, neighbors[1:5,])[, 1:3])
  
  expect_equal(neighbor_sum(annot$totalcounts, neighbors)[1:5],
               neighbor_sum(annot$totalcounts, neighbors[1:5, ]))
  
  expect_equal(neighbor_colSums(counts, neighbors)[1:5, ],
               neighbor_colSums(counts, neighbors[1:5, ]))
  
  expect_equal(neighbor_meantabulate(annot$celltype, neighbors)[1:5, 1:3],
               neighbor_meantabulate(annot$celltype, neighbors[1:5,])[, 1:3])
  
  expect_equal(neighbor_mean(annot$totalcounts, neighbors)[1:5],
               neighbor_mean(annot$totalcounts, neighbors[1:5, ]))
  
  expect_equal(neighbor_colMeans(counts, neighbors)[1:5, ],
               neighbor_colMeans(counts, neighbors[1:5, ]))
 
  expect_equal(neighbor_colMeans(counts, neighbors)[1:5, 1],
               neighbor_mean(counts[, 1], neighbors)[1:5])
 
  expect_equal(neighbor_colSums(counts, neighbors)[1:5, 1],
               neighbor_sum(counts[, 1], neighbors)[1:5])
  
  expect_equal(neighbor_colSums(counts, neighbors)[1:5, 1],
               neighbor_colSums(as.matrix(counts), neighbors)[1:5, 1])
})


#### test functions calling neighborhood summaries: ----------------
tmp <- dataframe_neighborhood_summary(df = annot[, c("fov", "totalcounts", "celltype")],
                                      neighbors = neighbors)
test_that("dataframe_neighborhood_summary returns a list of correct results", {
  expect_true(is.list(tmp))
  expect_true(is.vector(tmp$fov))
  expect_true(is.vector(tmp$totalcounts))
  expect_true(is.matrix(tmp$celltype))
  expect_true(dim(tmp$celltype)[2] == length(unique(annot$celltype)))
})

tmp <- get_neighborhood_expression(counts = counts, neighbors = neighbors)
test_that("get_neighborhood_expression is correct", {
  expect_identical(dim(tmp), dim(counts))
  expect_equal(unname(tmp[2, "SPP1"]),
               sum(counts[neighbors[2, ] > 0, "SPP1"]) / sum(neighbors[2, ] != 0))
})

#### test subsetting function: ----------------------------

test_that("neighbor subsampling works", {
  neighbors <- nearestNeighborGraph(x = cosmx_kidney$xy[, 1],
                                    y = cosmx_kidney$xy[, 2],
                                    N = 80)
  subsetted_neighbors <- subsampleNeighborsByRow(neighbors, p = 0.5)
  
  expect_true(all(Matrix::rowSums(subsetted_neighbors > 0) == 40)) # same number in each row
  expect_true(all(neighbors[subsetted_neighbors > 0] > 0))  # no new non-zero entries in subset
})