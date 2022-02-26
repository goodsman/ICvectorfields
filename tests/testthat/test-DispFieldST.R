test_that("DispFieldST correctly estimates vertical velocity", {
  Vec1 <- c(1:5, 0, 0, 0, 0)
  Mat1 <- Vec1
  for (i in 2:9) Mat1 <- rbind(Mat1, Vec1)

  Vec2 <- c(0, 1:5, 0, 0, 0)
  Mat2 <- Vec2
  for (i in 2:9) Mat2 <- rbind(Mat2, Vec2)

  Vec3 <- c(0, 0, 1:5, 0, 0)
  Mat3 <- Vec3
  for (i in 2:9) Mat3 <- rbind(Mat3, Vec3)

  Vec4 <- c(0, 0, 0, 1:5, 0)
  Mat4 <- Vec4
  for (i in 2:9) Mat4 <- rbind(Mat4, Vec4)

  rast1 <- terra::rast(Mat1)
  rast2 <- terra::rast(Mat2)
  rast3 <- terra::rast(Mat3)
  rast4 <- terra::rast(Mat4)

  teststack1 <- c(rast1, rast2, rast3, rast4)
  VFdf2 <- DispFieldST(teststack1, lag1 = 1, factv1 = 9, facth1 = 9)
  expect_equal(round(VFdf2$dispx, 7), 1)
  expect_equal(VFdf2$dispy, 0)
})
