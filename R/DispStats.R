#' Calculate statistics in source or sink regions
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
#' @param statistic desired output statistic: It should be one of "mean", "var",
#'   or "sum". Default setting is var.
#' @export
#'
#' @examples
#' # Illustrating use of DispStats:
#'
#' (Mat1 <- matrix(c(1,1,1,0,0,0,0,0,0,
#'                   1,1,1,0,0,0,0,0,0,
#'                   1,1,1,0,0,0,0,0,0,
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
#'                   1,1,1,0,0,0,0,0,0,
#'                   1,1,1,0,0,0,0,0,0,
#'                   1,1,1,0,0,0,0,0,0),
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
#' # Now to compute the statistics at the source: the mean of the original values
#' # in each region of interest (should be one in first row)
#' (VFdf2 <- DispStats(rast1, rast2, rast1, VFdf1, sourceloc = TRUE, statistic = "mean"))
#' # sum in each region of interest (should be nine in first row)
#' (VFdf3 <- DispStats(rast1, rast2, rast1, VFdf1, sourceloc = TRUE, statistic = "sum"))
#' # variance in each region of interest (should be zero in all rows)
#' (VFdf4 <- DispStats(rast1, rast2, rast1, VFdf1, sourceloc = TRUE, statistic = "var"))
DispStats <- function(inputrast1, inputrast2, statrast, vfdf,
                      sourceloc = TRUE, statistic = "var") {
  if (is.logical(sourceloc) == FALSE) {
    stop("sourceloc must be either TRUE or FALSE")
  }
  if (is.element(statistic, c("mean", "var", "sum")) == FALSE) {
    stop("statistic must be 'mean', 'var', or 'sum'")
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
      if (statistic == "mean") {
        Outdf$Mean[i] <- mean(inputmat3*mat2back*mat1bin, na.rm = TRUE)
      }
      if (statistic == "var") {
        Outdf$Var[i] <- stats::var(as.numeric(inputmat3*mat2back*mat1bin), na.rm = TRUE)
      }
      if (statistic == "sum") {
        Outdf$Sum[i] <- sum(inputmat3*mat2back*mat1bin, na.rm = TRUE)
      }
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
      if (statistic == "mean") {
        Outdf$Mean[i] <- mean(inputmat3*mat1forw*mat2bin, na.rm = TRUE)
      }
      if (statistic == "var") {
        Outdf$Var[i] <- stats::var(as.numeric(inputmat3*mat1forw*mat2bin), na.rm = TRUE)
      }
      if (statistic == "sum") {
        Outdf$Sum[i] <- sum(inputmat3*mat1forw*mat2bin, na.rm = TRUE)
      }
    }
  }
  return(Outdf)
}
