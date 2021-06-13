test_that("DispFieldST correctly estimates vertical velocity", {
  # creating matrices
  rseq <- stats::runif(54)
  Mat1 <- matrix(rep(0, 9 * 9), nrow = 9)
  Mat2 <- Mat1
  Mat3 <- Mat1
  Mat4 <- Mat1
  Mat1[1:9, 1:6] <- rseq
  Mat2[1:9, 2:7] <- rseq
  Mat3[1:9, 3:8] <- rseq
  Mat4[1:9, 4:9] <- rseq

  # Note that rasterizing a matrix causes it to be rotated 90 degrees.
  # Therefore, any shift in the x direction is in fact now a shift in the y direction
  rast1 <- terra::rast(Mat1)
  rast2 <- terra::rast(Mat2)
  rast3 <- terra::rast(Mat3)
  rast4 <- terra::rast(Mat4)

  teststack1 <- c(rast1, rast2, rast3, rast4)
  VFdf3 <- DispFieldSTbb(teststack1, lag1 = 1, 1, 9, 1, 9)
  expect_equal(round(VFdf3$dispy, 7), -0.1111111)
  expect_equal(VFdf3$dispx, 0)
})
