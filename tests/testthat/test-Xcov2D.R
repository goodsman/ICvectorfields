test_that("Xcov2D correctly returns matrix shifts", {
  Out1 <- ICvectorfields::GetRowCol(
    which.max(Xcov2D(matrix(c(1:6, rep(0, 3)), nrow = 3),
                     matrix(c(rep(0, 3), 1:6), nrow = 3))),
    dim1 = dim(Xcov2D(matrix(c(1:6, rep(0, 3)), nrow = 3),
                      matrix(c(rep(0, 3), 1:6), nrow = 3)))[1],
    dim2 = dim(Xcov2D(matrix(c(1:6, rep(0, 3)), nrow = 3),
                      matrix(c(rep(0, 3), 1:6), nrow = 3)))[2]
  )
  # See the examples section of Xcov2D function documentation for
  # an explanation of why c(6, 7) is expected.
  expect_equal(Out1, c(6, 7))
})
