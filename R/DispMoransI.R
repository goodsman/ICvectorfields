#' Calculate statistics in source or sink regions
#'
#' Functions for computing the statistics which may be driving variables of
#' movement that has been quantified using the \code{\link{DispField}} or
#' \code{\link{DispFieldbb}} functions. The same raster data as were supplied to
#' the aforementioned functions must be supplied to these in addition to a
#' raster layer for which statistics are sought. Then for each region of
#' interest defined when \code{\link{DispField}} or \code{\link{DispFieldbb}}
#' were called, these functions compute statistics for presumed source (sourceloc =
#' TRUE) locations or presumed sink locations (sourceloc = FALSE). Note that in the
#' DispMornasI function, defining radius using distance means that a radius of
#' one corresponds to the rook's neighbourhood.
#'
#' @rdname DispStats
#'
#' @param inputrast1 a raster as produced by terra::rast
#' @param inputrast2 a raster of equivalent dimension to inputrast1 as produced
#'   by terra::rast
#' @param statrast a raster of equivalent dimension to inputrast1 as produced by
#'   terra::rast which contains the variable that will be used to compute
#'   statistics
#' @param vfdf a data frame returned by the \code{\link{DispField}} or
#'   \code{\link{DispFieldbb}} functions, which contains all of the information
#'   necessary for defining regions of interest as well as the displacement
#'   estimates
#' @param sourceloc logical (TRUE or FALSE) indicating whether statistics are to be
#'   returned at source or sink locations
#' @param rad1 an ingeger indicating the neighbourhood radius for Moran's I
#'   statistic calculations in rows/columns. Any cell within a distance of rad1
#'   cells of the focal cell is considered to be in its neighbourhood.
#'
#' @return A data frame is returned with all of the same columns as the vfdf
#'   input data frame plus an additional column containing the computed
#'   statistic in each region of interest defined in vfdf.
#' @export
#'
#' @examples
#' # Illustrating use of DispMoransI:
#'
#' (Mat1 <- matrix(c(0.1,1,0.1,0,0,0,0,0,0,
#'                   1,0.1,1,0,0,0,0,0,0,
#'                   0.1,1,0.1,0,0,0,0,0,0,
#'                   0,0,0,0,0,0,0,0,0,
#'                   0,0,0,0,0,0,0,0,0,
#'                   0,0,0,0,0,0,0,0,0,
#'                   0,0,0,0,0,0,0,0,0,
#'                   0,0,0,0,0,0,0,0,0,
#'                   0,0,0,0,0,0,0,0,0),
#'                  nrow = 9))
#' (Mat2 <- matrix(c(0,0,0,0,0,0,0,0,0,
#'                   0,0,0,0,0,0,0,0,0,
#'                   0,0,0,0,0,0,0,0,0,
#'                   0,0,0,0,0,0,0,0,0,
#'                   0,0,0,0,0,0,0,0,0,
#'                   0,0,0,0,0,0,0,0,0,
#'                   0.1,1,0.1,0,0,0,0,0,0,
#'                   1,0.1,1,0,0,0,0,0,0,
#'                   0.1,1,0.1,0,0,0,0,0,0),
#'                  nrow = 9))
#'
#' # Note that rasterizing a matrix causes it to be rotated 90 degrees.
#' # Therefore, any shift in the x direction is in fact now a shift in the y direction
#' rast1 <- terra::rast(Mat1)
#' terra::plot(rast1)
#' rast2 <- terra::rast(Mat2)
#' terra::plot(rast2)
#'
#' (VFdf1 <- DispField(rast1, rast2, factv1 = 3, facth1 = 3))
#' # The second raster is shifted down by -0.6666667 units relative to the first raster
#' # dispy = -0.6666667 (the width of each box is 0.1111111).
#'
#' # Now to compute the statistics at the source: the Moran's I of the original values
#' # in each region of interest (should be minus one in first row)
#' (VFdf2 <- DispMoransI(rast1, rast2, rast1, VFdf1, sourceloc = TRUE, rad1 = 1))
#'
DispMoransI <- function(inputrast1, inputrast2, statrast, vfdf,
                      sourceloc = TRUE, rad1) {
  if (is.logical(sourceloc) == FALSE) {
    stop("sourceloc must be either TRUE or FALSE")
  }

  # Converting to matrix form
  inputmat1 <- RastToMatrix(inputrast1)
  inputmat2 <- RastToMatrix(inputrast2)
  inputmat3 <- RastToMatrix(statrast)

  # resolution variables
  dx <- terra::xres(inputrast1)
  dy <- terra::yres(inputrast1)

  Outdf <- vfdf

  # calculating shift in units of rows/columns
  shiftx = round(Outdf$dispx/dx)
  shifty = round(Outdf$dispy/dy)

  # preventing NA, Inf, and NAN shifts
  shiftx[is.na(shiftx) == TRUE] = 0.0
  shiftx[is.infinite(shiftx) == TRUE] = 0.0
  shiftx[is.nan(shiftx) == TRUE] = 0.0

  shifty[is.na(shifty) == TRUE] = 0.0
  shifty[is.infinite(shifty) == TRUE] = 0.0
  shifty[is.nan(shifty) == TRUE] = 0.0

  # cycling through all grid locations
  for (i in 1:dim(Outdf)[1]) {
    # extracting the relevant subset.
    mat1sub <- ExtractMat(inputmat1, Outdf$frowmin[i], Outdf$frowmax[i], Outdf$fcolmin[i], Outdf$fcolmax[i])
    mat2sub <- inputmat2

    # shifting and computing local statistics
    if (sourceloc == TRUE) {
      # converting source matrix to binary
      mat1bin <- mat1sub
      mat1bin[mat1bin > 0] <- 1
      mat1bin[mat1bin == 0] <- NA

      # shifting displaced matrix back
      if (abs(shifty[i]) < (dim(mat2sub)[1] - 1) &
          abs(shiftx[i]) < dim(mat2sub)[2] - 1) {
        mat2back <- ShiftMat(mat2sub,
                             shiftrows = -shifty[i],
                             shiftcols = -shiftx[i])
      } else {
        mat2back <- matrix(rep(0, dim(mat2sub)[1]*dim(mat2sub)[2]),
                           nrow = dim(mat2sub)[1])
      }

      # converting displaced matrix to binary
      mat2back[mat2back > 0] <- 1
      mat2back[mat2back == 0] <- NA

      # calculating the statistic
      # the na.rm argument is critical
      # Note that the variable for which
      # the statistic is sought is in inputmat3
      prodmat = inputmat3*mat2back*mat1bin
      MI <- MoransI(mat1 = prodmat, r1 = rad1)
      if (length(MI) == 0 || is.na(MI) || is.nan(MI) || is.infinite(MI) || MI == -999.0) {
        MI = NA
      }
      Outdf$MoransI[i] <- MI

    } else {
      # converting the second matrix to binary
      mat2bin <- mat2sub
      mat2bin[mat2bin > 0] <- 1
      mat2bin[mat2bin == 0] <- NA

      # shifting source matrix forward
      if (abs(shifty[i]) < (dim(mat2sub)[1] - 1) &
          abs(shiftx[i]) < dim(mat2sub)[2] - 1) {
        mat1forw <- ShiftMat(mat1sub,
                             shiftrows = shifty[i],
                             shiftcols = shiftx[i])
      } else {
        mat1forw <- matrix(rep(0, dim(mat2sub)[1]*dim(mat2sub)[2]),
                           nrow = dim(mat2sub)[1])
      }

      # converting displaced matrix to binary
      mat1forw[mat1forw > 0] <- 1
      mat1forw[mat1forw == 0] <- NA

      # calculating the statistic
      # the na.rm argument is critical
      # Note that the variable for which
      # the statistic is sought is in inputmat3
      prodmat = inputmat3*mat1forw*mat2bin
      MI <- MoransI(mat1 = prodmat, r1 = rad1)
      if (length(MI) == 0 || is.na(MI) || is.nan(MI) || is.infinite(MI) || MI == -999.0) {
        MI = NA
      }
    }
  }
  return(Outdf)
}
