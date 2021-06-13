test_that("FlipMat flips a matrix vertically and horizontally", {
  matz <- matrix(c(9, 6, 3, 8, 5, 2, 7, 4, 1),
                nrow = 3, byrow = T)
  expect_equal(FlipMat(matrix(c(1:9), nrow = 3)), matz)
})
