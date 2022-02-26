#' Diplacement fields using bounding box when velocity varies spatially
#'
#' This is an implementation of a novel algorithm that differs from more
#' traditional digital image correlation implementations that are applied in the
#' \code{\link{DispField}} and \code{\link{DispFieldbb}} functions. This version
#' is similar to the \code{\link{DispFieldSTbb}} function except that it does
#' not require a specific time lag. Instead the user specifies a maximum time
#' lag and the function computes displacement vectors using the time lag that
#' produces the maximum speed (magnitude of displacement divided by time lag).
#' The function calculates a displacement field representing persistent movement
#' based on the cross-covariance in a raster stack (in this case a sequential
#' series of rasters) presumably representing spatial population abundance or
#' density at more than two different instances of time. If analysis is
#' restricted to only two time instances, \code{\link{DispFieldbb}} is more
#' appropriate.
#'
#' The DispFieldSTbball function has the same inner workings as the
#' \code{\link{DispFieldSTbb}} function except that instead of specifying a
#' specific time lag, the user specifies a maximum time lag. The function then
#' cycles through all lags up to the maximum time lag and chooses the for each
#' location the maximum speed. The DispFieldSTbball function is more appropriate
#' than \code{\link{DispFieldSTbb}} when velocity is variable in space.
#'
#' Caution is warranted when defining the bounding box dimensions because the
#' function can produce erroneous results when the bounding box is too small.
#'
#' @param inputstack1 a raster stack with each raster layer representing an
#'   instance of time. The raster stack should be organized such that the first
#'   raster in the stack is the first observed spatial dataset and time
#'   progresses forward with the third dimension index of the raster stack. The
#'   raster stack should contain only numeric values. Any NA value will be
#'   converted to a zero
#' @param lagmax an integer representing the maximum time lag
#' @param rowmn an integer denoting the minimum row index of the sub-grid
#' @param rowmx an integer denoting the maximum row index of the sub-grid
#' @param colmn an integer denoting the minimum column index of the sub-grid
#' @param colmx an integer denoting the maximum column index of the sub-grid
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
#' @seealso \code{\link{DispFieldbb}} for a similar function with focal region
#'   defined using a bounding box for only two time instances,
#'   \code{\link{DispFieldSTbb}} for a version designed to quantify persistent
#'   directional movement when velocity is constant in space and the focal
#'   region is defined using a bounding box, see \code{\link{DispFieldSTall}}
#'   for a version designed to quantify persistent directional movement when
#'   velocity is variable in space and focal regions are defined based on a
#'   grid, and \code{\link{Xcov2D}} for demonstration of how two-dimensional
#'   cross-covariance is used to determine displacement (see examples of Xcov2D
#'   function documentation).
#'
#' @examples
#' (Mat1 <- matrix(rep(c(1:5, 0, 0, 0, 0), 9), nrow = 9, byrow = TRUE))
#' (Mat2 <- matrix(rep(c(0, 1:5, 0, 0, 0), 9), nrow = 9, byrow = TRUE))
#' (Mat3 <- matrix(rep(c(0, 0, 1:5, 0, 0), 9), nrow = 9, byrow = TRUE))
#' (Mat4 <- matrix(rep(c(0, 0, 0, 1:5, 0), 9), nrow = 9, byrow = TRUE))
#'
#' # Rasterizing
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
#' (VFdf5 <- DispFieldSTbball(teststack1, lagmax = 2, 1, 9, 1, 9))
#' # block is moving rightward at a speed of 1 unit of space per unit of time
#' # dispx = 1
#'
DispFieldSTbball <- function(inputstack1, lagmax, rowmn, rowmx, colmn, colmx, restricted = FALSE) {
  if (lagmax < 1 || lagmax != round(lagmax)) {
    stop("lagmax must be an integer larger than zero")
  }
  if (lagmax >= (dim(inputstack1)[3] - 1)) {
    stop("lagmax must be at least one smaller than the time demension")
  }
  if (lagmax == 1) {
    stop("no lag > 1; use DispFieldSTbb function")
  }
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

  # Initializing
  Outdf <- DispFieldSTbb(inputstack1, lag1 = 1, rowmn, rowmx, colmn, colmx, restricted = FALSE)

  # cycling through all larger lags
  for (i in 2:lagmax) {
    OutdfNew <- DispFieldSTbb(inputstack1, lag1 = i, rowmn, rowmx, colmn, colmx, restricted = FALSE)

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
