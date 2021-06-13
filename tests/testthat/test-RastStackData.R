test_that("RastStackData produces a raster stack with appropriate dimensions", {
  xyzdf <- expand.grid(x = c(1:3), y = c(1:3))
  xyzdf$z1 <- runif(9)
  xyzdf$z2 <- runif(9)
  xyzdf$z3 <- runif(9)

  zstack <- RastStackData(xyzdf)
  expect_equal(length(zstack), 3)
  expect_equal(dim(zstack[[1]]), c(3, 3, 1))
})
