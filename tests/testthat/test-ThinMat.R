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
