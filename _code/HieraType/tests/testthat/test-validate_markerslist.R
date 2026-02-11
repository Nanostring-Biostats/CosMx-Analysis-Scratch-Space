# Tests for validate_markerslist function

test_that("validate_markerslist removes missing predictors", {
  ml <- make_markerslist(
    index_marker = list(ct1 = "GENE1"),
    predictors = list(ct1 = c("GENE1", "GENE2", "GENE3", "GENE4", "GENE5"))
  )

  available <- c("GENE1", "GENE2", "GENE3")

  expect_warning(
    expect_message(
      result <- validate_markerslist(ml, available),
      "Removed.*predictor genes"
    ),
     "ct1 has only 3 positive predictor genes.  Consider adding more predictor genes for this celltype." 
  )

  
  expect_equal(sort(result$ct1$predictors), c("GENE1", "GENE2", "GENE3"))
  
})

test_that("validate_markerslist errors when all index markers missing", {
  ml <- make_markerslist(
    index_marker = list(ct1 = "GENE1"),
    predictors = list(ct1 = c("GENE1", "GENE2", "GENE3", "GENE4", "GENE5"))
  )

  available <- c("GENE2", "GENE3", "GENE4", "GENE5")

  expect_error(
    validate_markerslist(ml, available),
    "All specified index markers.*were not found"
  )
})

test_that("validate_markerslist warns when some index markers missing", {
  ml <- make_markerslist(
    index_marker = list(ct1 = c("GENE1", "GENE_MISSING")),
    predictors = list(ct1 = c("GENE1", "GENE2", "GENE3", "GENE4", "GENE5"))
  )

  available <- c("GENE1", "GENE2", "GENE3", "GENE4", "GENE5")

  expect_warning(
    result <- validate_markerslist(ml, available),
    "index markers.*were not found"
  )

  expect_equal(result$ct1$index_marker, "GENE1")
})

test_that("validate_markerslist warns about low predictor count", {
  ml <- make_markerslist(
    index_marker = list(ct1 = "GENE1"),
    predictors = list(ct1 = c("GENE1", "GENE2", "GENE3", "GENE4", "GENE5"))
  )

  available <- c("GENE1", "GENE2")

  expect_warning(
    validate_markerslist(ml, available),
    "has only.*positive predictor genes"
  )
})
