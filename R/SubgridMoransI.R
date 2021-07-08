#' Compute Moran's I for sub-grids
#'
#' A function that facilitates calculation of Moran's I at the sub-grid level.
#' This may be useful for evaluating spatial autocorrelation as a driver  of
#' movement speed or direction if used in tandem with \code{\link{DispField}},
#' \code{\link{DispFieldST}}, or \code{\link{DispFieldSTall}}.
#'
#' @importFrom Rcpp sourceCpp
#'
#' @param inputrast1 a raster as produced by terra::rast
#' @param factv1 an odd integer for the vertical dimension of sub-grids
#' @param facth1 an odd integer for the horizontal dimension of sub-grids
#' @param rad1 an integer indicating the neighbourhood radius in rows/columns.
#'   Any cell within a distance of rad1 cells of the focal cell is considered to
#'   be in its neighbourhood.
#'
#' @return A data frame is returned with the following column names: rowcent,
#'   colcent, frowmin, frowmax, fcolmin, fcolmax, and a column for the output
#'   Moran's I statistic.
#' @export
#'
#' @examples
#' (TestMat <- matrix(c(1, 0, 1, 0, 1,
#'                      0, 1, 0, 1, 0,
#'                      1, 0, 1, 0, 1,
#'                      0, 1, 0, 1, 0,
#'                      1, 0, 1, 0, 1),
#'                     nrow = 5))
#'
#' YestRast <- terra::rast(TestMat)
#' terra::plot(TestRast)
#'
#' SubgridMoransI(TestRast, factv1 = 5, facth1 = 5, rad1 = 1)
#' # using rad = 1 is equivalent to using the rooks neighbourhood
#' # and so the output should be -1.
#'
SubgridMoransI = function(inputrast1, factv1, facth1, rad1 = 1){
  if(factv1/2 == round(factv1/2)) {
    stop("factv1 and facth1 must be odd integers")
  }
  if(facth1/2 == round(facth1/2)) {
    stop("factv1 and facth1 must be odd integers")
  }
  if(round(factv1) != factv1 || round(facth1) != facth1) {
    stop("factv1 and facth1 must be integers")
  }
  if(rad1 > factv1 || rad1 > facth1) {
    stop("rad1 must be less than or equal to min(factv1, facth1)")
  }

  # converting output raster to matrix
  inputmat1 <- RastToMatrix(inputrast1)

  # obtaining the row and column indices for subsamples
  Outdf <- ThinMat(inputmat1, factv1, facth1)
  if(dim(Outdf)[1] < 1) stop("no viable grid locations: try smaller values
                             for factv1 and facth1")

  Outdf$MoransI <- rep(NA, dim(Outdf)[1])
  # cycling through all grid locations
  for(i in 1:dim(Outdf)[1]){
    mat1sub <- inputmat1[c(Outdf$frowmin[i]:Outdf$frowmax[i]), c(Outdf$fcolmin[i]:Outdf$fcolmax[i])]
    Outdf$MoransI[i] <- MoransI(mat1 = mat1sub, r1 = rad1)
  }

  return(Outdf)
}
