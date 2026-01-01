context("get.index function tests")

test_that("get.index correctly identifies root, nodes, and tips", {
  # Simple binary tree:
  #     3
  #    / \
  #   1   2
  edge <- matrix(c(3, 1,
                   3, 2), ncol = 2, byrow = TRUE)

  result <- get.index(edge)

  expect_equal(result$root, 3)
  expect_equal(sort(result$tips), c(1, 2))
  expect_equal(result$node, numeric(0))
})

test_that("get.index handles deeper trees", {
  # Tree structure:
  #       5
  #      / \
  #     4   3
  #    / \
  #   1   2
  edge <- matrix(c(5, 4,
                   5, 3,
                   4, 1,
                   4, 2), ncol = 2, byrow = TRUE)

  result <- get.index(edge)

  expect_equal(result$root, 5)
  expect_equal(sort(result$tips), c(1, 2, 3))
  expect_equal(result$node, 4)
})

test_that("get.index handles single edge tree", {
  edge <- matrix(c(2, 1), ncol = 2, byrow = TRUE)

  result <- get.index(edge)

  expect_equal(result$root, 2)
  expect_equal(result$tips, 1)
})
