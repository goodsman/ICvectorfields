#' Displacement fields based on 2D cross-covariance
#'
#' Calculates a displacement field based on the cross-covariance of two input
#' rasters presumably representing spatial population abundance or density at
#' two different instances of time.
#'
#' The input rasters are first converted to equivalent matrices. The function
#' then divides the domain up into sub-grids of size factv1 X facth1, which are
#' vertical and horizontal sub-grid dimensions.
#'
#' If restricted is set to FALSE (the default), the function computes
#' cross-covariance between each sub-grid of the first input raster and the
#' entirety of the second input raster and then uses the location of maximum
#' cross-covariance to estimate displacement in the vertical and horizontal
#' directions from the centre of each sub-grid.
#'
#' If restricted is set to TRUE, the function uses cross-covariance between each
#' sub-grid in the first input raster and the equivalent sub-grid in the second
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
#' Caution is warranted when defining the sub-grid dimensions because the
#' function can produce erroneous results when sub-grids are too small.
#'
#' @param inputrast1 a raster as produced by terra::rast
#' @param inputrast2 a raster of equivalent dimension to inputrast1 as produced
#'   by terra::rast
#' @param factv1 an odd integer for the vertical dimension of sub-grids
#' @param facth1 an odd integer for the horizontal dimension of sub-grids
#' @param restricted logical (TRUE or FALSE)
#'
#' @return A data frame is returned with the following column names: rowcent,
#'   colcent, frowmin, frowmax, fcolmin, fcolmax, centx, centy, dispx, and
#'   dispy. The rowcent and colcent column names are the row and column indices
#'   for the center of each sub-grid; frowmin and frowmax are the sub-grid
#'   minimum and maximum row indices; fcolmin and fcolmax are the sub-grid
#'   minimum and maximum column indices; centx and centy are the projected
#'   coordinates of the centre of the subgrid derived from the raster input
#'   files; dispx and dispy are the displacement in the horizontal and vertical
#'   directions in the same units as the projected coordinates of the raster
#'   input files.
#' @export
#'
#' @seealso \code{\link{DispFieldbb}} for a similar function using a bounding
#'   box to define a focal region, \code{\link{DispFieldST}} for a version
#'   designed to quantify persistent directional movement when the time series
#'   features more than two time instances, \code{\link{DispFieldSTall}} for a
#'   version designed to quantify persistent directional movement when velocity
#'   is variable in space, and \code{\link{Xcov2D}} for demonstration of how
#'   two-dimensional cross-covariance is used to determine displacement (see
#'   examples of Xcov2D function documentation).
#'
#' @examples
#' (Mat1 <- matrix(rep(c(1:5, 0, 0, 0, 0), 9), nrow = 9, byrow = TRUE))
#' (Mat2 <- matrix(rep(c(0, 1:5, 0, 0, 0), 9), nrow = 9, byrow = TRUE))
#'
#' # Note that rasterizing a matrix causes it to be rotated 90 degrees.
#' # Therefore, any shift in the x direction is in fact now a shift in the y direction
#' rast1 <- terra::rast(Mat1)
#' terra::plot(rast1)
#' rast2 <- terra::rast(Mat2)
#' terra::plot(rast2)
#'
#' (VFdf1 <- DispField(rast1, rast2, factv1 = 9, facth1 = 9))
#' # The second raster is shifted down by 0.1111111 units relative to the first raster
#' # dispy = -0.1111111
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
  dx <- terra::xres(inputrast1)
  dy <- terra::yres(inputrast1)

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
