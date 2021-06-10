test_that("ReflMat creates a relection matrix of dimension dim1 X dim1", {
  zmat = matrix(rep(0,9), nrow = 3, ncol = 3)
  zmat[1,3] = 1
  zmat[2,2] = 1
  zmat[3,1] = 1
  expect_equal(ReflMat(3), zmat)
})
