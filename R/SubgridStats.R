#' Compute statistics for subgrids
#'
#' @param inputrast1 a raster as produced by terra::rast
#' @param factv1 an odd integer for the vertical dimension of sub-grids
#' @param facth1 an odd integer for the horizontal dimension of sub-grids
#' @param statistic desired output statistic: It should be one of "mean",
#'   "var", or "sum". Default setting is var.
#' @export
#'
#' @examples
#' (TestMat <- matrix(c(1, 0, 1, 0, 1,
#'                     0, 1, 0, 1, 0,
#'                     1, 0, 1, 0, 1,
#'                     0, 1, 0, 1, 0,
#'                     1, 0, 1, 0, 1),
#'                     nrow = 5))
#'
#' TestRast <- terra::rast(TestMat)
#' terra::plot(TestRast)
#'
#' SubgridStats(TestRast, factv1 = 5, facth1 = 5, statistic = "mean")
#' SubgridStats(TestRast, factv1 = 5, facth1 = 5, statistic = "var")
#' SubgridStats(TestRast, factv1 = 5, facth1 = 5, statistic = "sum")
#'
SubgridStats <- function(inputrast1, factv1, facth1, statistic = "var") {
  if (factv1 / 2 == round(factv1 / 2)) {
    stop("factv1 and facth1 must be odd integers")
  }
  if (facth1 / 2 == round(facth1 / 2)) {
    stop("factv1 and facth1 must be odd integers")
  }
  if (round(factv1) != factv1 || round(facth1) != facth1) {
    stop("factv1 and facth1 must be integers")
  }
  if (is.element(statistic, c("mean", "var", "sum")) == FALSE) {
    stop("statistic must be 'mean', 'var', or 'sum'")
  }

  # converting output raster to matrix
  inputmat1 <- RastToMatrix(inputrast1)

  # obtaining the row and column indices for subsamples
  Outdf <- ThinMat(inputmat1, factv1, facth1)
  if (dim(Outdf)[1] < 1) stop("no viable grid locations: try smaller values
                             for factv1 and facth1")

  # adding a column for the statistic and computing it
  if (statistic == "mean") {
    Outdf$Mean <- rep(NA, dim(Outdf)[1])
    # cycling through all grid locations
    for (i in 1:dim(Outdf)[1]) {
      mat1sub <- inputmat1[c(Outdf$frowmin[i]:Outdf$frowmax[i]), c(Outdf$fcolmin[i]:Outdf$fcolmax[i])]
      Outdf$Mean[i] <- mean(mat1sub, na.rm = TRUE)
    }
  }
  if (statistic == "var") {
    Outdf$Var <- rep(NA, dim(Outdf)[1])
    # cycling through all grid locations
    for (i in 1:dim(Outdf)[1]) {
      mat1sub <- inputmat1[c(Outdf$frowmin[i]:Outdf$frowmax[i]), c(Outdf$fcolmin[i]:Outdf$fcolmax[i])]
      Outdf$Var[i] <- stats::var(as.numeric(mat1sub), na.rm = TRUE)
    }
  }
  if (statistic == "sum") {
    Outdf$Sum <- rep(NA, dim(Outdf)[1])
    for (i in 1:dim(Outdf)[1]) {
      mat1sub <- inputmat1[c(Outdf$frowmin[i]:Outdf$frowmax[i]), c(Outdf$fcolmin[i]:Outdf$fcolmax[i])]
      Outdf$Sum[i] <- sum(mat1sub, na.rm = TRUE)
    }
  }

  return(Outdf)
}
