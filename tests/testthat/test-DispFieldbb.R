test_that("DispFieldbb correctly estimates displacement", {
  rseq <- stats::runif(72)
  Mat1 <- matrix(rep(0, 81), nrow = 9)
  Mat2 <- Mat1
  Mat1[1:9, 1:8] <- rseq
  Mat2[1:9, 2:9] <- rseq

  rast1 <- terra::rast(Mat1)
  rast2 <- terra::rast(Mat2)

  VFdf1 <- DispFieldbb(rast1, rast2, 1, 9, 1, 9)
  expect_equal(round(VFdf1$dispx, 7), 0.1111111)
  expect_equal(VFdf1$dispy, 0)
})
