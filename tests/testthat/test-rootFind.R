context("rootFind function tests")

test_that("rootFind correctly identifies root node", {
  # Simple tree with root = 3
  edge <- matrix(c(3, 1,
                   3, 2), ncol = 2, byrow = TRUE)

  expect_equal(rootFind(edge), 3)
})

test_that("rootFind works with deeper trees", {
  # Tree with root = 5
  edge <- matrix(c(5, 4,
                   5, 3,
                   4, 1,
                   4, 2), ncol = 2, byrow = TRUE)

  expect_equal(rootFind(edge), 5)
})

test_that("rootFind errors when no unique root exists", {
  # Cyclic structure - no unique root
  edge <- matrix(c(1, 2,
                   2, 3,
                   3, 1), ncol = 2, byrow = TRUE)

  expect_error(rootFind(edge))
})
