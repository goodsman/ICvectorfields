test_that("RotationDetect correctly classifies rotation", {
  # creating rotation patterns
  Mat1 <- matrix(rep(0,9*9), nrow = 9)
  Mat1[c(1:3), 4] <- 1
  Mat1[c(7:9), 6] <- 1
  Mat1[4, c(7:9)] <- 1
  Mat1[6, c(1:3)] <- 1

  Mat2 <- matrix(rep(0,9*9), nrow = 9)
  Mat2[c(1:3), 5] <- 1
  Mat2[c(7:9), 5] <- 1
  Mat2[5, c(7:9)] <- 1
  Mat2[5, c(1:3)] <- 1

  # rasterizing
  rast1 <- terra::rast(Mat1)
  rast2 <- terra::rast(Mat2)

  # Detecting counter-clockwise rotation
  VFdf1 <- DispField(rast1, rast2, factv1 = 3, facth1 = 3, restricted = TRUE)
  patdf1 <- RotationDetect(VFdf1)
  expect_equal(patdf1$Pattern[5], "clockwise")

  # Detecting clockwise rotation
  VFdf2 <- DispField(rast2, rast1, factv1 = 3, facth1 = 3, restricted = TRUE)
  patdf2 <- RotationDetect(VFdf2)
  expect_equal(patdf2$Pattern[5], "counter")
})
