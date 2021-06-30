test_that("PixelCt function correctly counts non-zero pixels", {
  Mat1 <- matrix(rep(c(1:5, 0, 0, 0, 0), 9), nrow = 9, byrow = TRUE)
  Mat2 <- matrix(rep(c(0, 1:5, 0, 0, 0), 9), nrow = 9, byrow = TRUE)

  # Note that rasterizing a matrix causes it to be rotated 90 degrees.
  # Therefore, any shift in the x direction is in fact now a shift in the y direction
  rast1 <- terra::rast(Mat1)
  rast2 <- terra::rast(Mat2)

  Confdf1 <- PixelCt(c(rast1, rast2), factv1 = 9, facth1 = 9)
  expect_equal(Confdf1$pixelct[1], 54)
})
