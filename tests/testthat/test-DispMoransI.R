test_that("DispMoransI correctly computes MoransI at the source", {
  Mat1 <- matrix(c(0.1,1,0.1,0,0,0,0,0,0,
                    1,0.1,1,0,0,0,0,0,0,
                    0.1,1,0.1,0,0,0,0,0,0,
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
                    0.1,1,0.1,0,0,0,0,0,0,
                    1,0.1,1,0,0,0,0,0,0,
                    0.1,1,0.1,0,0,0,0,0,0),
                  nrow = 9)

  # Note that rasterizing a matrix causes it to be rotated 90 degrees.
  # Therefore, any shift in the x direction is in fact now a shift in the y direction
  rast1 <- terra::rast(Mat1)
  rast2 <- terra::rast(Mat2)

  VFdf1 <- DispField(rast1, rast2, factv1 = 3, facth1 = 3)
  # The second raster is shifted down by -0.6666667 units relative to the first raster
  # dispy = -0.6666667 (the width of each box is 0.1111111).

  # Now to compute the statistics at the source: the Moran's I of the original values
  # in each region of interest (should be minus one in first row)
  VFdf2 <- DispMoransI(rast1, rast2, rast1, VFdf1, sourceloc = TRUE, rad1 = 1)
  expect_equal(round(VFdf2$MoransI[1], 3), -1)
})
