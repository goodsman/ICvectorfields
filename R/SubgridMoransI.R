#' Compute Moran's I for sub-grids
#'
#' @importFrom Rcpp sourceCpp
#'
#' @param inputrast1
#' @param factv1
#' @param facth1
#' @param rad1
#'
#' @return
#' @export
#'
#' @examples
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
    Outdf$MoransI[i] <- MoransI(mat1 = mat1sub, rad1 = rad1)
  }

  return(Outdf)
}
