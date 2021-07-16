test_that("DispStats computes shifted statistics correctly", {
  Mat1 <- matrix(c(1,1,1,0,0,0,0,0,0,
                   1,1,1,0,0,0,0,0,0,
                   1,1,1,0,0,0,0,0,0,
                   0,0,0,0,0,0,0,0,0,
                   0,0,0,0,0,0,0,0,0,
                   0,0,0,0,0,0,0,0,0,
                   0,0,0,0,0,0,0,0,0,
                   0,0,0,0,0,0,0,0,0,
                   0,0,0,0,0,0,0,0,0),
                  nrow = 9)
  Mat2 <- matrix(c(0,0,0,0,0,0,0,0,0,
                   0,0,0,0,0,0,0,0,0,
                   0,0,0,0,0,0,0,0,0,
                   0,0,0,0,0,0,0,0,0,
                   0,0,0,0,0,0,0,0,0,
                   0,0,0,0,0,0,0,0,0,
                   1,1,1,0,0,0,0,0,0,
                   1,1,1,0,0,0,0,0,0,
                   1,1,1,0,0,0,0,0,0),
                  nrow = 9)

  # Note that rasterizing a matrix causes it to be rotated 90 degrees.
  # Therefore, any shift in the x direction is in fact now a shift in the y direction
  rast1 <- terra::rast(Mat1)
  rast2 <- terra::rast(Mat2)

  # Calculating displacement
  VFdf1 <- DispField(rast1, rast2, factv1 = 3, facth1 = 3)

  # Calculating statistics at the source
  VFdf2 <- DispStats(rast1, rast2, rast1, VFdf1, sourceloc = TRUE, statistic = "mean")
  VFdf3 <- DispStats(rast1, rast2, rast1, VFdf1, sourceloc = TRUE, statistic = "sum")
  VFdf4 <- DispStats(rast1, rast2, rast1, VFdf1, sourceloc = TRUE, statistic = "var")

  expect_equal(VFdf2$Mean[1], 1)
  expect_equal(VFdf3$Sum[1], 9)
  expect_equal(VFdf4$Var[1], 0)
})
