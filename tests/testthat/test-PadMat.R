test_that("PadMat returns a square matrix with even number of rows and coluns", {
  matz <- matrix(rep(0, 100), nrow = 10)
  expect_equal(PadMat(matrix(rep(0, 6), nrow = 2)), matz)
})
