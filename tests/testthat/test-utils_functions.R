test_that("ThinMat works", {
  rc <- 2
  cc <- 2
  frmn <- 1
  frmx <- 3
  fcmn <- 1
  fcmx <- 3
  dfz <- data.frame(rc, cc, frmn, frmx, fcmn, fcmx)
  colnames(dfz) <- c(
    "rowcent", "colcent", "frowmin", "frowmax",
    "fcolmin", "fcolmax"
  )
  expect_equal(
    ThinMat(matrix(c(1:9), nrow = 3, ncol = 3), factv = 3, facth = 3),
    dfz
  )
})

test_that("ReflMat function in utils.R outputs correct matrix", {
  zmat <- matrix(rep(0, 9), nrow = 3)
  zmat[1, 3] <- 1
  zmat[2, 2] <- 1
  zmat[3, 1] <- 1
  expect_equal(ReflMat(3), zmat)
})

test_that("raster to matrix conversion works", {
  rastz <- terra::rast(matrix(1:9, nrow = 3))
  matz <- matrix(1:9, nrow = 3, byrow = TRUE)
  expect_equal(RastToMatrix(rastz), matz)
})

test_that("PadMat returns a square matrix with even number of rows and coluns", {
  matz <- matrix(rep(0, 100), nrow = 10)
  expect_equal(PadMat(matrix(rep(0, 6), nrow = 2)), matz)
})

test_that("FlipMat flips a matrix vertically and horizontally", {
  matz <- matrix(c(9, 6, 3, 8, 5, 2, 7, 4, 1),
                 nrow = 3, byrow = T
  )
  expect_equal(FlipMat(matrix(c(1:9), nrow = 3)), matz)
})

test_that("ExtractMat extracts relevant elements and sets the rest to zero", {
  matz <- matrix(c(rep(0, 5), 6, 10, 0, 0, 7, 11, rep(0, 5)), nrow = 4, byrow = TRUE)
  expect_equal(ExtractMat(matrix(c(1:16), nrow = 4), rowmin = 2, rowmax = 3, colmin = 2, colmax = 3), matz)
})
