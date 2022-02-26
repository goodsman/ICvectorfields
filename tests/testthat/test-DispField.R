test_that("DispField correctly estimates displacement", {
  # constructing test matrices
  Vec1 <- c(1:5, 0, 0, 0, 0)
  Mat1 <- Vec1
  for (i in 2:9) Mat1 <- rbind(Mat1, Vec1)

  Vec2 <- c(0, 1:5, 0, 0, 0)
  Mat2 <- Vec2
  for (i in 2:9) Mat2 <- rbind(Mat2, Vec2)

  # converting to rasters
  rast1 <- terra::rast(Mat1)
  rast2 <- terra::rast(Mat2)

  VFdf1 <- DispField(rast1, rast2, factv1 = 9, facth1 = 9)
  expect_equal(round(VFdf1$dispx, 7), 1)
  expect_equal(VFdf1$dispy, 0)
})
