# Tests for make_markerslist function

test_that("make_markerslist creates valid markerslist object", {
  ml <- make_markerslist(
    index_marker = list(ct1 = "GENE1", ct2 = "GENE2"),
    predictors = list(ct1 = c("GENE1", "GENE3", "GENE4"),
                      ct2 = c("GENE2", "GENE5", "GENE6"))
  )

  expect_s3_class(ml, "markerslist")
  expect_equal(length(ml), 2)
  expect_equal(names(ml), c("ct1", "ct2"))
})

test_that("make_markerslist validates missing names in index_marker", {
  expect_error(
    make_markerslist(
      index_marker = list(c("GENE1")),
      predictors = list(ct1 = c("GENE1", "GENE2"))
    ),
    "index_marker must be a named list"
  )
})

test_that("make_markerslist validates missing names in predictors", {
  expect_error(
    make_markerslist(
      index_marker = list(ct1 = "GENE1"),
      predictors = list(c("GENE1", "GENE2"))
    ),
    "predictors must be a named list"
  )
})

test_that("make_markerslist validates mismatched celltype names", {
  expect_error(
    make_markerslist(
      index_marker = list(ct1 = "GENE1"),
      predictors = list(ct2 = c("GENE1", "GENE2"))
    ),
    "Some of the celltype class names"
  )
})

test_that("make_markerslist adds index_marker to predictors if missing", {
  expect_message(
    ml <- make_markerslist(
      index_marker = list(ct1 = "GENE1"),
      predictors = list(ct1 = c("GENE2", "GENE3"))
    ),
    "will be added to the predictors list"
  )

  expect_true("GENE1" %in% ml$ct1$predictors)
})

test_that("make_markerslist handles use_offclass_markers_as_negative_predictors", {
  ml <- make_markerslist(
    index_marker = list(ct1 = "GENE1", ct2 = "GENE2"),
    predictors = list(ct1 = c("GENE1", "GENE3"),
                      ct2 = c("GENE2", "GENE4")),
    use_offclass_markers_as_negative_predictors = FALSE
  )

  expect_false(ml$ct1$use_offclass_markers_as_negative_predictors)
  expect_false(ml$ct2$use_offclass_markers_as_negative_predictors)
})
