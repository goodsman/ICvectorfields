#' Displacement fields based on 2D cross-covariance using bounding box
#'
#' Calculates a displacement field based on the cross-covariance of two input
#' rasters presumably representing spatial population abundance or density at
#' two different instances of time. This version differs from
#' \code{\link{DispField}} in that the user defines a bounding box that
#' determines a single sub-grid. The center of the bounding box is the location
#' from whence displacement is estimated.
#'
#' The input rasters are first converted to equivalent matrices. If restricted
#' is set to FALSE (the default), the function computes cross-covariance between
#' the sub-grid of the first input raster and the entirety of the second input
#' raster and then uses the location of maximum cross-covariance to estimate
#' displacement in the vertical and horizontal directions from the centre of the
#' sub-grid.
#'
#' If restricted is set to TRUE, the function uses cross-covariance between the
#' sub-grid of the first input raster and the equivalent sub-grid of the second
#' input raster to estimate vertical and horizontal displacement.
#'
#' Reference coordinates and cell size are extracted from the first input raster
#' such that the locations from whence displacement is estimated as well as
#' displacement estimates can be expressed in the units of the projected
#' coordinates.
#'
#' The coordinates are assumed to increase vertically and horizontally from the
#' lower left corner of the two-dimensional domain.
#'
#' Caution is warranted when defining the bounding box because the function can
#' produce erroneous results when the bounding box is too small.
#'
#' @param inputrast1 a raster as produced by terra::rast
#' @param inputrast2 a raster of equivalent dimension to inputrast1 as produced
#'   by terra::rast
#' @param rowmn an integer denoting the minimum row index of the sub-grid
#' @param rowmx an integer denoting the maximum row index of the sub-grid
#' @param colmn an integer denoting the minimum column index of the sub-grid
#' @param colmx an integer denoting the maximum column index of the sub-grid
#' @param restricted logical (TRUE or FALSE)
#'
#' @return A data frame is returned with the following column names: rowcent,
#'   colcent, frowmin, frowmax, fcolmin, fcolmax, centx, centy, dispx, and
#'   dispy. The rowcent and colcent column names are the row and column indices
#'   for the center of the sub-grid; frowmin and frowmax are the sub-grid
#'   minimum and maximum row indices; fcolmin and fcolmax are the sub-grid
#'   minimum and maximum column indices; centx and centy are the projected
#'   coordinates of the centre of the subgrid derived from the raster input
#'   files; dispx and dispy are the displacement in the horizontal and vertical
#'   directions in the same units as the projected coordinates of the raster
#'   input files.
#' @export
#'
#' @seealso \code{\link{DispField}} for a similar function with a grid of focal
#'   regions, \code{\link{DispFieldSTbb}} for a version designed to quantify
#'   persistent directional movement when the time series features more than two
#'   time instances, \code{\link{DispFieldSTbball}} for a version designed to
#'   quantify persistent directional movement when velocity is variable in
#'   space, and \code{\link{Xcov2D}} for demonstration of how two-dimensional
#'   cross-covariance is used to determine displacement (see examples of Xcov2D
#'   function documentation).
#'
#' @examples
#' rseq <- stats::runif(72)
#' Mat1 <- matrix(rep(0, 81), nrow = 9)
#' Mat2 <- Mat1
#' Mat1[1:9, 1:8] <- rseq
#' Mat1
#' Mat2[1:9, 2:9] <- rseq
#' Mat2
#'
#' # rasterizing
#' rast1 <- terra::rast(Mat1)
#' terra::plot(rast1)
#' rast2 <- terra::rast(Mat2)
#' terra::plot(rast2)
#'
#' (VFdf1 <- DispFieldbb(rast1, rast2, 2, 8, 2, 8))
#' # The second raster is shifted right by 1 unit relative to the first raster
#' # dispx = 1
DispFieldbb <- function(inputrast1, inputrast2, rowmn, rowmx, colmn, colmx, restricted = FALSE) {
  if (rowmn != round(rowmn)) {
    stop("rowmn and rowmx must be positive integers")
  }
  if (rowmx != round(rowmx)) {
    stop("rowmn and rowmx must be positive integers")
  }
  if (colmn != round(colmn)) {
    stop("colmn and colmx must be positive integers")
  }
  if (colmx != round(colmx)) {
    stop("colmn and colmx must be positive integers")
  }
  if (is.logical(restricted) == FALSE) {
    stop("restricted must be either TRUE or FALSE")
  }

  # Converting to matrix form
  inputmat1 <- RastToMatrix(inputrast1)
  inputmat2 <- RastToMatrix(inputrast2)

  # Obtaining the row and column indices for subsample
  rowcent <- rowmn + floor((rowmx - rowmn) / 2)
  colcent <- colmn + floor((colmx - colmn) / 2)
  Outdf <- data.frame(rowcent, colcent)
  Outdf$frowmin <- rowmn
  Outdf$frowmax <- rowmx
  Outdf$fcolmin <- colmn
  Outdf$fcolmax <- colmx

  # Adding columns for central coordinates and displacement
  Outdf$centx <- rep(NA, dim(Outdf)[1])
  Outdf$centy <- rep(NA, dim(Outdf)[1])
  Outdf$dispx <- rep(NA, dim(Outdf)[1])
  Outdf$dispy <- rep(NA, dim(Outdf)[1])

  # resolution variables
  dx <- terra::xres(inputrast1)
  dy <- terra::yres(inputrast1)

  # cycling through all grid locations
  for (i in 1:dim(Outdf)[1]) {
    if (restricted == FALSE) {
      # extracting the relevant subset.
      mat1sub <- ExtractMat(inputmat1, rowmn, rowmx, colmn, colmx)
      mat2sub <- inputmat2
    } else {
      mat1sub <- inputmat1[c(rowmn:rowmx), c(colmn:colmx)]
      mat2sub <- inputmat2[c(rowmn:rowmx), c(colmn:colmx)]
    }

    # computing cross-covariance
    if (sum(mat1sub) > 0) {
      # XcovMat = Xcov2D(mat1sub, inputmat2)
      XcovMat <- Xcov2D(mat1sub, mat2sub)
      # Computing displacement vector
      xcoord1 <- GetRowCol(which.max(XcovMat), dim1 = dim(XcovMat)[1], dim2 = dim(XcovMat)[2])[2]
      ycoord1 <- GetRowCol(which.max(XcovMat), dim1 = dim(XcovMat)[1], dim2 = dim(XcovMat)[2])[1]
      # translate rows and columns to coordinates
      Outdf$centx[i] <- terra::xFromCol(inputrast1, col = Outdf$colcent[i])
      Outdf$centy[i] <- terra::yFromRow(inputrast1, row = Outdf$rowcent[i])
      # Computing displacement: because this is a square matrix with an even number of rows,
      # the center is at the center.
      Outdf$dispx[i] <- (xcoord1 - (dim(XcovMat)[1] / 2 + 1)) * dx
      Outdf$dispy[i] <- ((dim(XcovMat)[2] / 2 + 1) - ycoord1) * dy
    } else {
      Outdf$centx[i] <- terra::xFromCol(inputrast1, col = Outdf$colcent[i])
      Outdf$centy[i] <- terra::yFromRow(inputrast1, row = Outdf$rowcent[i])
      Outdf$dispx[i] <- 0
      Outdf$dispy[i] <- 0
    }
  }
  return(Outdf)
}
