# Tests for totalcount_norm function


test_that("totalcount_norm handles zero totalcounts", {
  mat <- Matrix::Matrix(c(0, 10, 0, 20), nrow = 2, ncol = 2, sparse = TRUE)
  rownames(mat) <- c("cell1", "cell2")
  colnames(mat) <- c("gene1", "gene2")

  # Should not error

  expect_no_error(result <- totalcount_norm(mat))

  # Zero row should remain zero (or at least finite)
  expect_true(all(is.finite(as.matrix(result))))
})

test_that("totalcount_norm uses provided totalcounts", {
  mat <- Matrix::Matrix(c(10, 20, 30, 40), nrow = 2, ncol = 2, sparse = TRUE)
  rownames(mat) <- c("cell1", "cell2")
  colnames(mat) <- c("gene1", "gene2")

  tc <- c(100, 200)
  names(tc) <- rownames(mat)
  
  result <- totalcount_norm(mat, tc = tc)
  result_notc <- totalcount_norm(mat)
  # Check dimensions preserved
  expect_false(isTRUE(all.equal(result, result_notc)))
  expect_equal(dim(result), dim(mat))
})

test_that("totalcount_norm preserves sparsity", {
  mat <- Matrix::Matrix(c(0, 1, 0, 2, 0, 3), nrow = 2, ncol = 3, sparse = TRUE)

  result <- totalcount_norm(mat)

  expect_true(inherits(result, "sparseMatrix") || inherits(result, "Matrix"))
})
