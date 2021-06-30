#' Count populated pixels in a raster stack
#'
#' @param inputstack1
#' @param factv1
#' @param facth1
#'
#' @return
#' @export
#'
#' @examples
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
  Outdf <- ThinMat(inputmat1, factv1, facth1)
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
