test_that("SubgridMoransI computes Moran's I correctly", {
  TestMat <- matrix(c(1, 0, 1, 0, 1,
                      0, 1, 0, 1, 0,
                      1, 0, 1, 0, 1,
                      0, 1, 0, 1, 0,
                      1, 0, 1, 0, 1),
                     nrow = 5)

  TestRast <- terra::rast(TestMat)
  Testdf <- SubgridMoransI(TestRast, factv1 = 5, facth1 = 5, rad1 = 1)
  expect_equal(round(Testdf$MoransI[1], 3), -1)
})
