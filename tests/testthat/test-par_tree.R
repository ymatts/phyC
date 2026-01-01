context("par.tree function tests")

test_that("par.tree validates input data type", {
  # Create a data frame with non-numeric columns
  vaf <- data.frame(
    Normal = c("a", "b", "c"),
    Sample1 = c("d", "e", "f")
  )

  expect_error(
    par.tree(vaf),
    "Invalid data type"
  )
})

test_that("par.tree validates VAF values range", {
  # Create VAF data with values outside valid range
  vaf <- data.frame(
    Normal = c(0.1, 0.2, 0.3),
    Sample1 = c(1.5, 0.5, 0.3)  # 1.5 is outside [0, 1]
  )

  expect_error(
    par.tree(vaf),
    "VAF values must be between 0 and 1"
  )
})

test_that("par.tree validates minimum columns", {
  # Single column - need at least Normal + 1 sample
  vaf <- data.frame(Normal = c(0.1, 0.2, 0.3))

  expect_error(
    par.tree(vaf),
    "VAF data must have at least 2 columns"
  )
})

test_that("par.tree handles valid input correctly", {
  skip_on_cran()

  # Create valid VAF data
  vaf <- data.frame(
    Normal = c(0, 0, 0, 0),
    Sample1 = c(0.3, 0.4, 0.1, 0.2),
    Sample2 = c(0.2, 0.3, 0.4, 0.1)
  )

  result <- par.tree(vaf)

  expect_true(is.list(result))
  expect_true("tree" %in% names(result))
  expect_true("edge" %in% names(result))
  expect_true("edge.length" %in% names(result))
})
