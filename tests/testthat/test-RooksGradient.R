test_that("RooksGradient correctly computes the gradient", {
  # creating pattern patterns
  Mat1 <- matrix(rep(0,9*9), nrow = 9)
  Mat1[c(4:6), c(4:6)] <- 2
  Mat1[c(4:6), c(1:3)] <- 1
  Mat1[c(1:3), c(4:6)] <- 1
  Mat1[c(7:9), c(4:6)] <- 1
  Mat1[c(4:6), c(7:9)] <- 1

  Rast1 <- terra::rast(Mat1)

  # calculating the mean in 9 subgrids
  statsdf1 <- SubgridStats(Rast1, factv1 = 3, facth1 = 3, statistic = "mean")

  # computing the gradient statistic on the mean
  graddf1 <- RooksGradient(statsdf1, statistic = "mean")

  expect_equal(graddf1$Gradient[1], 1)
  expect_equal(graddf1$Gradient[5], -1)
})
