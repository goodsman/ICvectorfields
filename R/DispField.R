DispField <- function(inputrast1, inputrast2, factv1, facth1, restricted = FALSE) {
  if (factv1 / 2 == round(factv1 / 2)) {
    stop("factv1 and facth1 must be odd integers")
  }
  if (facth1 / 2 == round(facth1 / 2)) {
    stop("factv1 and facth1 must be odd integers")
  }
  if (round(factv1) != factv1 || round(facth1) != facth1) {
    stop("factv1 and facth1 must be integers")
  }
  if (is.logical(restricted) == FALSE) {
    stop("restricted must be either TRUE or FALSE")
  }

  # Converting to matrix form
  inputmat1 <- RastToMatrix(inputrast1)
  inputmat2 <- RastToMatrix(inputrast2)

  # Obtaining the row and column indices for subsamples
  Outdf <- ThinMat(inputmat1, factv1, facth1)
  if (dim(Outdf)[1] < 1) stop("no viable grid locations: try smaller values
                             for factv1 and facth1")

  # Adding columns for central coordinates and displacement
  Outdf$centx <- rep(NA, dim(Outdf)[1])
  Outdf$centy <- rep(NA, dim(Outdf)[1])
  Outdf$dispx <- rep(NA, dim(Outdf)[1])
  Outdf$dispy <- rep(NA, dim(Outdf)[1])

  # resolution variables
  dx <- xres(inputrast1)
  dy <- yres(inputrast1)

  # cycling through all grid locations
  for (i in 1:dim(Outdf)[1]) {
    if (restricted == FALSE) {
      # extracting the relevant subset.
      mat1sub <- ExtractMat(inputmat1, Outdf$frowmin[i], Outdf$frowmax[i], Outdf$fcolmin[i], Outdf$fcolmax[i])
      mat2sub <- inputmat2
    } else {
      mat1sub <- inputmat1[c(Outdf$frowmin[i]:Outdf$frowmax[i]), c(Outdf$fcolmin[i]:Outdf$fcolmax[i])]
      mat2sub <- inputmat2[c(Outdf$frowmin[i]:Outdf$frowmax[i]), c(Outdf$fcolmin[i]:Outdf$fcolmax[i])]
    }

    # computing cross-covariance
    if (sum(mat1sub) > 0) {
      # XcovMat = Xcov2D(mat1sub, inputmat2)
      XcovMat <- Xcov2D(mat1sub, mat2sub)
      # Computing displacement vector
      xcoord1 <- GetRowCol(which.max(XcovMat), dim1 = dim(XcovMat)[1], dim2 = dim(XcovMat)[2])[2]
      ycoord1 <- GetRowCol(which.max(XcovMat), dim1 = dim(XcovMat)[1], dim2 = dim(XcovMat)[2])[1]
      # translate rows and columns to coordinates
      Outdf$centx[i] <- xFromCol(inputrast1, col = Outdf$colcent[i])
      Outdf$centy[i] <- yFromRow(inputrast1, row = Outdf$rowcent[i])
      # Computing displacement: because this is a square matrix with an even number of rows,
      # the center is at the center.
      Outdf$dispx[i] <- (xcoord1 - (dim(XcovMat)[1] / 2 + 1)) * dx
      Outdf$dispy[i] <- ((dim(XcovMat)[2] / 2 + 1) - ycoord1) * dy
    } else {
      Outdf$centx[i] <- xFromCol(inputrast1, col = Outdf$colcent[i])
      Outdf$centy[i] <- yFromRow(inputrast1, row = Outdf$rowcent[i])
      Outdf$dispx[i] <- 0
      Outdf$dispy[i] <- 0
    }
  }
  return(Outdf)
}
