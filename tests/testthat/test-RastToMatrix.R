test_that("raster to matrix conversion works", {
  rastz <- terra::rast(matrix(1:9, nrow = 3))
  matz <- matrix(1:9, nrow = 3, byrow = TRUE)
  expect_equal(RastToMatrix(rastz), matz)
})
