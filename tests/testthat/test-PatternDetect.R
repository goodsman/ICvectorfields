test_that("PatternDetect correctly classifies patterns", {
  # matrices to create convergence/divergence patterns
  Mat1 <- matrix(rep(0, 9*9), nrow = 9)
  Mat1[3, c(4, 6)] <- 1
  Mat1[7, c(4, 6)] <- 1
  Mat1[c(4, 6), 3] <- 1
  Mat1[c(4, 6), 7] <- 1

  Mat2 <- matrix(rep(0, 9*9), nrow = 9)
  Mat2[2, c(4, 6)] <- 1
  Mat2[8, c(4, 6)] <- 1
  Mat2[c(4, 6), 2] <- 1
  Mat2[c(4, 6), 8] <- 1

  # rasterizing
  rast1 <- terra::rast(Mat1)
  rast2 <- terra::rast(Mat2)

  # Detecting a divergence
  VFdf1 <- DispField(rast1, rast2, factv1 = 3, facth1 = 3, restricted = TRUE)
  patdf1 <- PatternDetect(VFdf1)
  expect_equal(patdf1$Pattern[5], "divergence")

  # Detecting a convergence
  VFdf2 <- DispField(rast2, rast1, factv1 = 3, facth1 = 3, restricted = TRUE)
  patdf2 <- PatternDetect(VFdf2)
  expect_equal(patdf2$Pattern[5], "convergence")
})
