#' Compute statistics for subgrids
#'
#' Functions that facilitate calculation of statistics at the sub-grid level.
#' These may be useful for drivers  of movement speed or direction if used in
#' tandem with \code{\link{DispField}}, \code{\link{DispFieldST}}, or
#' \code{\link{DispFieldSTall}}.
#'
#' Note that when using radius to define the neighbourhood in Moran's I
#' calculations, a radius of one corresponds to the rook's neibhourhood. Values
#' that are NA or Inf are not included in calculations of the Moran's I
#' statistic nor in any of the other statistics that can be computed.
#'
#' @rdname SubgridStats
#' @param inputrast1 a raster as produced by terra::rast
#' @param factv1 an odd integer for the vertical dimension of sub-grids
#' @param facth1 an odd integer for the horizontal dimension of sub-grids
#' @param rad1 an integer indicating the neighbourhood radius for Moran's I
#'   statistic calculations in rows/columns. Any cell within a distance of rad1
#'   cells of the focal cell is considered to be in its neighbourhood.
#'
#' @return A data frame is returned with the following column names: rowcent,
#'   colcent, frowmin, frowmax, fcolmin, fcolmax, and a column for the output
#'   statistic.
#' @export
#'
#' @seealso \code{\link{DispStats}} and \code{\link{DispMoransI}}for functions
#'   that compute statistics at presumed source or sink locations in each region
#'   of interest.
#'
#' @examples
#' (TestMat <- matrix(c(1, 0, 1, 0, 1,
#'                      0, 1, 0, 1, 0,
#'                      1, 0, 1, 0, 1,
#'                      0, 1, 0, 1, 0,
#'                      1, 0, 1, 0, 1),
#'                     nrow = 5))
#'
#' TestRast <- terra::rast(TestMat)
#' terra::plot(TestRast)
#'
#' SubgridMoransI(TestRast, factv1 = 5, facth1 = 5, rad1 = 1)
#' # using rad1 = 1 is equivalent to using the rooks neighbourhood
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
  # note it is very important to not simply
  # convert NA/Inf values to zeros when calculating
  # Moran's I because this greatly effects the statistic.
  # Instead, NA/Inf values should be flagged as NA
  # so that the MoransI function can disregard them.
  inputmat1 <- RastToMatrix(inputrast1, NAapproach = "NA")

  # obtaining the row and column indices for subsamples
  Outdf <- ThinMat(inputmat1, factv1, facth1)
  if(dim(Outdf)[1] < 1) stop("no viable grid locations: try smaller values
                             for factv1 and facth1")

  Outdf$MoransI <- rep(NA, dim(Outdf)[1])
  # cycling through all grid locations
  for(i in 1:dim(Outdf)[1]){
    mat1sub <- inputmat1[c(Outdf$frowmin[i]:Outdf$frowmax[i]), c(Outdf$fcolmin[i]:Outdf$fcolmax[i])]
    MI <- MoransI(mat1 = mat1sub, r1 = rad1)
    if (length(MI) == 0 || is.na(MI) || is.nan(MI) || is.infinite(MI) || MI == -999.0) {
      MI = NA
    }
    Outdf$MoransI[i] <- MI
  }

  return(Outdf)
}
