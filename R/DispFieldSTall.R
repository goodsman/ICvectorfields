#' Diplacement fields for spatiotemporal data when velocity varies spatially
#'
#' This is an implementation of a novel algorithm that differs from more
#' traditional digital image correlation implementations that are applied in the
#' \code{\link{DispField}} and \code{\link{DispFieldbb}} functions. This version
#' is similar to the \code{\link{DispFieldST}} function except that it does not
#' require a specific time lag. Instead the user specifies a maximum time lag
#' and the function computes displacement vectors using the time lag that
#' produces the maximum speed (magnitude of displacement divided by time lag).
#' The function calculates a displacement field representing persistent movement
#' based on the cross-covariance in a raster stack (in this case a sequential
#' series of rasters) presumably representing spatial population abundance or
#' density at more than two different instances of time. If analysis is
#' restricted to only two time instances, \code{\link{DispField}} is more
#' appropriate.
#'
#' The DispFieldSTall function has the same inner workings as the
#' \code{\link{DispFieldST}} function except that instead of specifying a
#' specific time lag, the user specifies a maximum time lag. The function then
#' cycles through all lags up to the maximum time lag and choses the for each
#' location the maximum speed. The DispFieldSTall function is more appropriate
#' than \code{\link{DispFieldST}} when velocity is variable in space.
#'
#' Caution is warranted when defining the sub-grid dimensions because the
#' function can produce erroneous results when sub-grids are too small.
#'
#' @param inputstack1 a raster stack with each raster layer representing an
#'   instance of time. The raster stack should be organized such that the first
#'   raster in the stack is the first observed spatial dataset and time
#'   progresses forward with the third dimension index of the raster stack. The
#'   raster stack should contain only numeric values. Any NA value will be
#'   converted to a zero
#' @param lagmax an integer representing the maximum time lag
#' @param factv1 an odd integer for the vertical dimension of subgrids
#' @param facth1 an odd integer for the horizontal dimension of subgrids
#' @param restricted logical (TRUE or FALSE)
#'
#' @return A data frame is returned with the following column names: rowcent,
#'   colcent, frowmin, frowmax, fcolmin, fcolmax, centx, centy, dispx, and
#'   dispy. The rowcent and colcent column names are the row and column indices
#'   for the center of each sub-grid; frowmin and frowmax are the sub-grid
#'   minimum and maximum row indices; fcolmin and fcolmax are the sub-grid
#'   minimum and maximum column indices; centx and centy are the projected
#'   coordinates of the centre of the subgrid derived from the raster input
#'   files; dispx and dispy are the orthoganal velocity vectors in units of
#'   space per timestep in the horizontal and vertical directions in the same
#'   spatial units as the projected coordinates of the raster input files.
#' @export
#'
#' @seealso \code{\link{DispField}} for a similar function with a grid of focal
#'   regions for only two time instances, \code{\link{DispFieldST}} for a
#'   version designed to quantify persistent directional movement when the time
#'   series features more than two time instances and the velocity is constant
#'   in space, \code{\link{DispFieldSTbball}} for a version designed to quantify
#'   persistent directional movement when velocity is variable in space and the
#'   focal region is defined using a bounding box, and \code{\link{Xcov2D}} for
#'   demonstration of how two-dimensional cross-covariance is used to determine
#'   displacement (see examples of Xcov2D function documentation).
#'
#' @examples
#' (Mat1 <- matrix(rep(c(1:5, 0, 0, 0, 0), 9), nrow = 9, byrow = TRUE))
#' (Mat2 <- matrix(rep(c(0, 1:5, 0, 0, 0), 9), nrow = 9, byrow = TRUE))
#' (Mat3 <- matrix(rep(c(0, 0, 1:5, 0, 0), 9), nrow = 9, byrow = TRUE))
#' (Mat4 <- matrix(rep(c(0, 0, 0, 1:5, 0), 9), nrow = 9, byrow = TRUE))
#'
#' # Note that rasterizing a matrix causes it to be rotated 90 degrees.
#' # Therefore, any shift in the x direction is in fact now a shift in the y direction
#' rast1 <- terra::rast(Mat1)
#' terra::plot(rast1)
#' rast2 <- terra::rast(Mat2)
#' terra::plot(rast2)
#' rast3 <- terra::rast(Mat3)
#' terra::plot(rast3)
#' rast4 <- terra::rast(Mat4)
#' terra::plot(rast4)
#'
#' teststack1 <- c(rast1, rast2, rast3, rast4)
#' (VFdf4 <- DispFieldSTall(teststack1, lagmax = 2, factv1 = 9, facth1 = 9))
#' # block is moving downward at a speed of 0.1111111 units of space per unit of time
#' # dispy = -0.1111111
DispFieldSTall <- function(inputstack1, lagmax, factv1, facth1, restricted = FALSE) {
  if (lagmax < 1 || lagmax != round(lagmax)) {
    stop("lagmax must be an integer larger than zero")
  }
  if (lagmax >= (dim(inputstack1)[3] - 1)) {
    stop("lagmax must be at least two smaller than the time demension")
  }
  if (lagmax == 1) {
    stop("no lag > 1; use DispFieldST function")
  }
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

  # Initializing
  Outdf <- DispFieldST(inputstack1, lag1 = 1, factv1, facth1, restricted = restricted)

  # cycling through all larger lags
  for (i in 2:lagmax) {
    OutdfNew <- DispFieldST(inputstack1, lag1 = i, factv1, facth1, restricted = restricted)

    # If the magnitude of displacement in OutdfNew > Outdf, then dispx and
    # dispy are replaced with the new values
    MagNew <- sqrt((OutdfNew$dispx^2) + OutdfNew$dispy^2)
    MagOld <- sqrt((Outdf$dispx^2) + Outdf$dispy^2)
    MagNew[is.na(MagNew) == TRUE] <- 0
    MagOld[is.na(MagOld) == TRUE] <- 0
    Outdf$dispx[MagNew > MagOld] <- OutdfNew$dispx[MagNew > MagOld]
    Outdf$dispy[MagNew > MagOld] <- OutdfNew$dispy[MagNew > MagOld]
  }

  return(Outdf)
}
