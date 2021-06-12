test_that("DispFieldbb correctly estimates displacement", {
  rseq = stats::runif(72)
  Mat1 = matrix(rep(0, 81), nrow = 9)
  Mat2 = Mat1
  Mat1[1:9, 1:8] = rseq
  Mat1
  Mat2[1:9, 2:9] = rseq
  Mat2

  # Note that rasterizing a matrix causes it to be rotated 90 degrees.
  # Therefore, any shift in the x direction is in fact now a shift in the y direction
  rast1 = terra::rast(Mat1)
  terra::plot(rast1)
  rast2 = terra::rast(Mat2)
  terra::plot(rast2)

  VFdf1 = DispFieldbb(rast1, rast2, 2, 8, 2, 8)
  expect_equal(round(VFdf1$dispy, 7), -0.1111111)
  expect_equal(VFdf1$dispx, 0)
})
