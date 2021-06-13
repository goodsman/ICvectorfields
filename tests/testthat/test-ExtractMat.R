test_that("ExtractMat extracts relevant elements and sets the rest to zero", {
  matz <- matrix(c(rep(0, 5), 6, 10, 0, 0, 7, 11, rep(0, 5)), nrow = 4, byrow = TRUE)
  expect_equal(ExtractMat(matrix(c(1:16), nrow = 4), rowmin = 2, rowmax = 3, colmin = 2, colmax = 3), matz)
})
