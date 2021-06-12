test_that("row and column retrival works", {
  # Assumes that the elements of the matrix are filled by column (byrow =
  # FALSE), which is the default
  expect_equal(GetRowCol(6, dim1 = 3, dim2 = 3), c(3, 2))
})
