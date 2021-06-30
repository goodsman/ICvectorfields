#' Count populated pixels in a raster stack
#'
#' In order to produce reliable results, the cross-covariance approach
#' implemented in \code{\link{DispField}}. \code{\link{DispFieldST}}, and
#' \code{\link{DispFieldSTall}} needs a certain minimum number of non-zero or
#' non-NA valued pixels in pairs of images or pairs of arrays derived from a
#' raster stack that it uses to compute cross-covariance. The user may define a
#' threshold such as ten percent of the pixels within each sub-grid. This
#' function allows the user to assess whether the minimum threshold number of
#' non-zero pixels per sub-grid are surpassed by returning the number of
#' non-zero pixels within each sub-grid over all of the time instances in the
#' user-supplied raster stack. The user may choose to disregard or mistrust
#' displacement or velocity estimates derived from sub-grids with insufficient
#' numbers of non-zero pixels. This function is designed to be called before or
#' after one of the functions referenced above in order to enable the user to
#' quantify their confidence in displacement or velocity estimates.
#'
#' @param inputstack1 a raster stack with each raster layer representing an
#'   instance of time. The raster stack should be organized such that the first
#'   raster in the stack is the first observed spatial dataset and time
#'   progresses forward with the third dimension index of the raster stack. The
#'   raster stack should contain only numeric values. Any NA value will be
#'   converted to a zero
#' @param factv1 an odd integer for the vertical dimension of subgrids
#' @param facth1 an odd integer for the horizontal dimension of subgrids
#'
#' @return A data frame is returned with the following column names: rowcent,
#'   colcent, frowmin, frowmax, fcolmin, fcolmax, and PixelCt. The rowcent and
#'   colcent column names are the row and column indices for the center of each
#'   sub-grid; frowmin and frowmax are the sub-grid minimum and maximum row
#'   indices; fcolmin and fcolmax are the sub-grid minimum and maximum column
#'   indices; pixelct is the count of non-zero pixels in the sub-grid over the
#'   entire time period covered by the input raster stack.
#' @export
#'
#' @examples
#' # below the example in the DispField function documentation is reproduced
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
#' (Confdf1 <- PixelCt(c(rast1, rast2), factv1 = 9, facth1 = 9))
#' # This should return a pixel count of 54: This is the number
#' # of pixels that were occupied in either the first or second
#' # time instance.
PixelCt <- function(inputstack1, factv1, facth1) {
  if (factv1 / 2 == round(factv1 / 2)) {
    stop("factv1 and facth1 must be odd integers")
  }
  if (facth1 / 2 == round(facth1 / 2)) {
    stop("factv1 and facth1 must be odd integers")
  }
  if (round(factv1) != factv1 || round(facth1) != facth1) {
    stop("factv1 and facth1 must be integers")
  }

  # Summing each layer and converting values greater than zero to one.
  SummedStack <- sum(inputstack1, na.rm = TRUE)
  SummedStack[SummedStack > 0] <- 1

  # converting output raster to matrix
  Summedmatrix <- RastToMatrix(SummedStack)

  # obtaining the row and column indices for subsamples
  Outdf <- ThinMat(Summedmatrix, factv1, facth1)
  if (dim(Outdf)[1] < 1) stop("no viable grid locations: try smaller values
                             for factv1 and facth1")

  # adding a column for pixel counts
  Outdf$pixelct <- rep(NA, dim(Outdf)[1])

  # cycling through all grid locations
  for (i in 1:dim(Outdf)[1]) {
    mat1sub <- Summedmatrix[c(Outdf$frowmin[i]:Outdf$frowmax[i]), c(Outdf$fcolmin[i]:Outdf$fcolmax[i])]
    Outdf$pixelct[i] <- sum(mat1sub)
  }

  return(Outdf)
}
