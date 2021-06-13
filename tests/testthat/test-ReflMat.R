test_that("ReflMat function in utils.R outputs correct matrix", {
  zmat <- matrix(rep(0, 9), nrow = 3)
  zmat[1, 3] <- 1
  zmat[2, 2] <- 1
  zmat[3, 1] <- 1
  expect_equal(ReflMat(3), zmat)
})
