#' Displacement fields for spatiotemporal data using a bounding box
#'
#' This is an implementation of a novel algorithm that differs from more
#' traditional digital image correlation implementations that are applied in the
#' \code{\link{DispField}} and \code{\link{DispFieldbb}} functions. The function
#' calculates a displacement field representing persistent movement based on the
#' cross-covariance in a raster stack (in this case a sequential series of
#' rasters) presumably representing spatial population abundance or density at
#' more than two different instances of time. If analysis is restricted to only
#' two time instances, \code{\link{DispFieldbb}} is more appropriate.
#'
#' The input rasters in the raster stack are first converted to equivalent
#' matrices, which together represent a three-dimensional array with two spatial
#' dimensions and one time dimension. The prescribed lag is applied to the three
#' dimensional array derived from the raster stack by first producing two
#' equivalent arrays and then removing appropriate numbers of layers from the
#' top of one and the bottom of the other. These are referred to as unlagged and
#' lagged spatiotemporal arrays in the description that follows.
#'
#' The function converts three dimensional lagged and unlagged spatiotemporal
#' arrays to two-dimensional lagged and unlagged spatiotemporal matrices by
#' averaging along one of the spatial dimensions (either rows or columns) to
#' obtain two pairs of two-dimensional matrices in which one dimension is
#' spatial (either rows or columns) and one dimension is temporal. One of each
#' pair corresponds to the  unlagged spatiotemporal array and the other
#' corresponds to the lagged spatiotemporal array. Displacement in the vertical
#' direction is computed using unlagged and lagged matrices that have been
#' averaged along rows and displacement in the horizontal direction is computed
#' using unlagged and lagged matrices that have been averaged along columns.
#'
#' If restricted is set to FALSE (the default), the function computes
#' cross-covariance between the values within the bounding box of the unlagged
#' row-averaged spatiotemporal matrix and the whole row-averaged lagged
#' spatiotemporal matrix and between the values within the bounding box of the
#' unlagged column-averaged spatiotemporal matrix and the entirety corresponding
#' lagged matrix.
#'
#' If restricted is set to TRUE, the function uses cross-covariance between
#' lagged and unlagged version of row-averaged and column averaged
#' spatiotemporal matrices that have all been either row or column-averaged
#' within the bounding box to estimate vertical and horizontal displacement.
#'
#' Regardless of whether restricted is set to TRUE or FALSE, for each sub-grid,
#' displacement in the x and y direction is divided by the shift in the time
#' dimension to produce orthogonal velocity vetors. Note that for this reason,
#' the lag1 argument of the function does not necessarily determine the time lag
#' that is used to produce each orthoganal velocity vector.
#'
#' Reference coordinates and cell size are extracted from the first raster stack
#' such that the locations from whence displacement is estimated as well as
#' displacement (or velocity) estimates can be expressed in the units of the
#' projected coordinates.
#'
#' The coordinates are assumed to increase vertically and horizontally from the
#' lower left corner of the two-dimensional domain.
#'
#' Caution is warranted when defining the sub-grid dimensions because the
#' function can produce erroneous results when sub-grids are too small.
#'
#' #' In addition, results can be quite sensitive to specification of the time
#' lag. If velocities are highly variable in space or over time, avoid
#' specifying a single time lag by calling the related
#' \code{\link{DispFieldSTbball}} function.
#'
#' @param inputstack1 a raster stack with each raster layer representing an
#'   instance of time. The raster stack should be organized such that the first
#'   raster in the stack is the first observed spatial dataset and time
#'   progresses forward with the third dimension index of the raster stack. The
#'   raster stack should contain only numeric values. Any NA value will be
#'   converted to a zero
#' @param lag1 an integer time lag
#' @param rowmn an integer for the minimum row index of the bounding box
#' @param rowmx an integer for the maximum row index of the bounding box
#' @param colmn an integer for the minimum column index of the bounding box
#' @param colmx an integer for the maximum column index of the bounding box
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
#'   series features more than two time instances but using a grid to define
#'   focal regions, see \code{\link{DispFieldSTbball}} for a version designed to
#'   quantify persistent directional movement when velocity is variable in
#'   space, and \code{\link{Xcov2D}} for demonstration of how two-dimensional
#'   cross-covariance is used to determine displacement (see examples of Xcov2D
#'   function documentation).
#'
#' @examples
#' rseq <- stats::runif(54)
#' Mat1 <- matrix(rep(0, 9*9), nrow = 9)
#' Mat2 <- Mat1; Mat3 <- Mat1; Mat4 <- Mat1
#' Mat1[1:9, 1:6] <- rseq
#' Mat1
#' Mat2[1:9, 2:7] <- rseq
#' Mat2
#' Mat3[1:9, 3:8] <- rseq
#' Mat3
#' Mat4[1:9, 4:9] <- rseq
#' Mat4
#'
#' # rasterizing
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
#' (VFdf3 <- DispFieldSTbb(teststack1, lag1 = 1, 2, 8, 2, 8))
#' # block is moving rightward at a speed of 0.1111111 units of space per unit of time
#' # dispx = 0.1111111
DispFieldSTbb <- function(inputstack1, lag1, rowmn, rowmx, colmn, colmx, restricted = FALSE) {
  if (lag1 < 1 || lag1 != round(lag1)) {
    stop("lag1 must be an integer larger than zero")
  }
  if (lag1 >= (dim(inputstack1)[3] - 1)) {
    stop("lag must be at least two smaller than the time demension")
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

  # breaking the stack into component rasters and converting to matrix form
  # and converting NA values to zero.
  for (i in 1:dim(inputstack1)[3]) assign(paste("inputmat", i, sep = ""), RastToMatrix(inputstack1[[i]]))

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
  dx <- terra::xres(inputstack1)
  dy <- terra::yres(inputstack1)

  # cycling through all grid locations
  for (i in 1:dim(Outdf)[1]) {
    if (restricted == FALSE) {
      # Constructing the space time matrices
      focalx <- c()
      focaly <- c()
      bufferx <- c()
      buffery <- c()
      for (j in 1:dim(inputstack1)[3]) {
        focalx <- rbind(focalx, colMeans(ExtractMat(
          get(paste("inputmat", j, sep = "")),
          rowmn, rowmx, colmn, colmx
        )))
        focaly <- rbind(focaly, rowMeans(ExtractMat(
          get(paste("inputmat", j, sep = "")),
          rowmn, rowmx, colmn, colmx
        )))
        bufferx <- rbind(bufferx, colMeans(get(paste("inputmat", j, sep = ""))))
        buffery <- rbind(buffery, rowMeans(get(paste("inputmat", j, sep = ""))))
      }

      # applying the lag
      lag1seq <- 1:lag1
      revlag1seq <- (dim(inputstack1)[3] - lag1) + lag1seq
      bufferx <- bufferx[-lag1seq, ]
      buffery <- buffery[-lag1seq, ]
      focalx <- focalx[-revlag1seq, ]
      focaly <- focaly[-revlag1seq, ]
    } else {
      # Constructing the space time matrices
      focalx <- c()
      focaly <- c()
      for (j in 1:dim(inputstack1)[3]) {
        focalx <- rbind(focalx, colMeans(get(paste("inputmat", j, sep = ""))[c(rowmn:rowmx), c(colmn:colmx)]))
        focaly <- rbind(focaly, rowMeans(get(paste("inputmat", j, sep = ""))[c(rowmn:rowmx), c(colmn:colmx)]))
      }

      # applying the lag
      lag1seq <- 1:lag1
      revlag1seq <- (dim(inputstack1)[3] - lag1) + lag1seq
      bufferx <- focalx[-lag1seq, ]
      buffery <- focaly[-lag1seq, ]
      focalx <- focalx[-revlag1seq, ]
      focaly <- focaly[-revlag1seq, ]
    }

    # computing cross-covariance and displacement in the x direction
    if (sum(focalx) > 0) {
      XcovMat <- Xcov2D(focalx, bufferx)
      # Computing displacement vector
      xcoord1 <- GetRowCol(which.max(XcovMat), dim1 = dim(XcovMat)[1], dim2 = dim(XcovMat)[2])[2]
      # the time coordinate is along margin 1 (rows or y)
      tcoord1 <- GetRowCol(which.max(XcovMat), dim1 = dim(XcovMat)[1], dim2 = dim(XcovMat)[2])[1]
      # translate rows and columns to coordinates
      Outdf$centx[i] <- terra::xFromCol(inputstack1[[1]], col = Outdf$colcent[i])
      Outdf$centy[i] <- terra::yFromRow(inputstack1[[1]], row = Outdf$rowcent[i])
      # Computing displacement: because this is a square matrix with an even number of rows,
      # the center is at the center.
      dispt <- ((dim(XcovMat)[2] / 2 + 1) - tcoord1) + lag1
      if (dispt != 0) Outdf$dispx[i] <- (xcoord1 - (dim(XcovMat)[1] / 2 + 1)) * dx / dispt
    } else {
      Outdf$centx[i] <- terra::xFromCol(inputstack1[[1]], col = Outdf$colcent[i])
      Outdf$centy[i] <- terra::yFromRow(inputstack1[[1]], row = Outdf$rowcent[i])
      Outdf$dispx[i] <- 0
    }

    # computing cross-covariance and displacement in the y direction
    if (sum(focaly) > 0) {
      XcovMat <- Xcov2D(focaly, buffery)
      # Computing displacement vector. Below I use the column index because space runs in
      # the horizontal direction. This is different from the non-spatiotemporal version.
      ycoord1 <- GetRowCol(which.max(XcovMat), dim1 = dim(XcovMat)[1], dim2 = dim(XcovMat)[2])[2]
      # the time coordinate is along margin 1 (rows or y)
      tcoord1 <- GetRowCol(which.max(XcovMat), dim1 = dim(XcovMat)[1], dim2 = dim(XcovMat)[2])[1]
      # translate rows and columns to coordinates
      Outdf$centx[i] <- terra::xFromCol(inputstack1[[1]], col = Outdf$colcent[i])
      Outdf$centy[i] <- terra::yFromRow(inputstack1[[1]], row = Outdf$rowcent[i])
      # Computing displacement: because this is a square matrix with an even number of rows,
      # the center is at the center.
      dispt <- ((dim(XcovMat)[2] / 2 + 1) - tcoord1) + lag1
      if (dispt != 0) Outdf$dispy[i] <- -(ycoord1 - (dim(XcovMat)[1] / 2 + 1)) * dy / dispt
    } else {
      Outdf$centx[i] <- terra::xFromCol(inputstack1[[1]], col = Outdf$colcent[i])
      Outdf$centy[i] <- terra::yFromRow(inputstack1[[1]], row = Outdf$rowcent[i])
      Outdf$dispy[i] <- 0
    }
  } # This ends the i loop

  return(Outdf)
}
