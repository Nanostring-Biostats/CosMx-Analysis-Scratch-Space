# Tests for make_pipeline function

test_that("make_pipeline creates valid pipeline object", {
  ml1 <- make_markerslist(
    index_marker = list(ct1 = "GENE1", ct2 = "GENE2"),
    predictors = list(ct1 = c("GENE1", "GENE3"), ct2 = c("GENE2", "GENE4"))
  )

  pipeline <- make_pipeline(markerslists = list("level1" = ml1))

  expect_s3_class(pipeline, "pipeline")
  expect_equal(names(pipeline), c("markerslists", "priors", "priors_category"))
})

test_that("make_pipeline validates that markerslists have names", {
  ml1 <- make_markerslist(
    index_marker = list(ct1 = "GENE1"),
    predictors = list(ct1 = c("GENE1", "GENE2"))
  )

  expect_error(
    make_pipeline(markerslists = list(ml1)),
    regexp = NULL
  )
})

test_that("make_pipeline validates priors category celltype reference existing markerslists", {
  ml1 <- make_markerslist(
    index_marker = list(ct1 = "GENE1", ct2 = "GENE2"),
    predictors = list(ct1 = c("GENE1", "GENE3"), ct2 = c("GENE2", "GENE4"))
  )
  ml2 <- make_markerslist(
    index_marker = list(ct1_sub = "GENE5"),
    predictors = list(ct1_sub = c("GENE5", "GENE6"))
  )

  expect_error(
    make_pipeline(
      markerslists = list("level1" = ml1, "level2" = ml2),
      priors = list("level2" = "nonexistent"),
      priors_category = list("level2" = "ct1_sub")
    ),
    "`priors_category` 'ct1_sub' not found as a celltype in the prior markerslist 'nonexistent'"
  )
  
})

test_that("make_pipeline validates priors category celltype reference existing markerslists, v2", {
  ml1 <- make_markerslist(
    index_marker = list(ct1 = "GENE1", ct2 = "GENE2"),
    predictors = list(ct1 = c("GENE1", "GENE3"), ct2 = c("GENE2", "GENE4"))
  )
  ml2 <- make_markerslist(
    index_marker = list(ct1_sub = "GENE5"),
    predictors = list(ct1_sub = c("GENE5", "GENE6"))
  )

  
  expect_error(
    make_pipeline(
      markerslists = list("level1" = ml1, "level2" = ml2),
      priors = list("level2" = "ct1"), ## should be "level1"
      priors_category = list("level2" = "ct1")
    ),
    "`priors_category` 'ct1' not found as a celltype in the prior markerslist 'ct1'"
    )
})

test_that("make_pipeline validates priors_category references valid celltypes", {
  ml1 <- make_markerslist(
    index_marker = list(ct1 = "GENE1", ct2 = "GENE2"),
    predictors = list(ct1 = c("GENE1", "GENE3"), ct2 = c("GENE2", "GENE4"))
  )
  ml2 <- make_markerslist(
    index_marker = list(ct1_sub = "GENE5"),
    predictors = list(ct1_sub = c("GENE5", "GENE6"))
  )

  expect_error(
    make_pipeline(
      markerslists = list("level1" = ml1, "level2" = ml2),
      priors = list("level2" = "level1"),
      priors_category = list("level2" = "nonexistent_celltype")
    ),
    "not found as a celltype in the prior markerslist"
  )
})

test_that("make_pipeline requires at least one root markerslist", {
  ml1 <- make_markerslist(
    index_marker = list(ct1 = "GENE1"),
    predictors = list(ct1 = c("GENE1", "GENE2"))
  )

  expect_error(
    make_pipeline(
      markerslists = list("level1" = ml1),
      priors = list("level1" = "level1"),
      priors_category = list("level1" = "ct1")
    ),
    "At least one celltype should be omitted from priors"
  )
})

test_that("make_pipeline accepts hierarchical structure", {
  ml1 <- make_markerslist(
    index_marker = list(ct1 = "GENE1", ct2 = "GENE2"),
    predictors = list(ct1 = c("GENE1", "GENE3"), ct2 = c("GENE2", "GENE4"))
  )
  ml2 <- make_markerslist(
    index_marker = list(ct1_sub1 = "GENE5", ct1_sub2 = "GENE6"),
    predictors = list(ct1_sub1 = c("GENE5", "GENE7"), ct1_sub2 = c("GENE6", "GENE8"))
  )

  pipeline <- make_pipeline(
    markerslists = list("level1" = ml1, "level2" = ml2),
    priors = list("level2" = "level1"),
    priors_category = list("level2" = "ct1")
  )

  expect_s3_class(pipeline, "pipeline")
  expect_equal(pipeline$priors$level2, "level1")
  expect_equal(pipeline$priors_category$level2, "ct1")
})
