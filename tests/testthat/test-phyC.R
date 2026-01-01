context("phyC main function tests")

test_that("phyC validates input parameters", {
  # Test with invalid cluster number
  edgeList <- list(
    matrix(c(3, 1, 3, 2), ncol = 2, byrow = TRUE)
  )
  edgeLenList <- list(c(1, 1))

  # Cluster number exceeds number of trees

  expect_error(
    phyC(edgeList, edgeLenList, cluster = 5),
    "cluster number cannot exceed"
  )
})

test_that("phyC validates edgeList format", {
  # Test with NULL input
  expect_error(
    phyC(NULL, list(c(1, 1))),
    "edgeList must be a non-empty list"
  )

  # Test with empty list
  expect_error(
    phyC(list(), list()),
    "edgeList must be a non-empty list"
  )
})

test_that("phyC validates edgeLenList format", {
  edgeList <- list(
    matrix(c(3, 1, 3, 2), ncol = 2, byrow = TRUE)
  )

  # Test with mismatched lengths
  expect_error(
    phyC(edgeList, list(c(1, 1), c(2, 2))),
    "edgeList and edgeLenList must have the same length"
  )
})

test_that("phyC handles valid input correctly", {
  skip_on_cran()

  # Create simple test data with at least 2 trees for clustering
  edgeList <- list(
    matrix(c(3, 1, 3, 2), ncol = 2, byrow = TRUE),
    matrix(c(3, 1, 3, 2), ncol = 2, byrow = TRUE)
  )
  edgeLenList <- list(c(1, 1), c(2, 2))

  result <- phyC(edgeList, edgeLenList, cluster = 2)

  expect_true(is.list(result))
  expect_true("trees" %in% names(result))
  expect_true("cluster" %in% names(result))
  expect_true("dist" %in% names(result))
})
