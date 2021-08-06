#' Calculate Gradient Statistics in the Rook's Neighbourhood
#'
#' The movement of populations into adjacent cells may sometimes be influenced
#' by the gradient of some predictive variable. This function enables the
#' calculation of a simple gradient statistic in the rook's neighbourhood of
#' each cell in a dataset. The statistic must first be computed for each grid
#' cell using \code{\link{SubgridStats}}. Then for each grid, the RooksGradient
#' function computes the arithmetic average of the difference between the
#' statistic at the focal grid cell and the statistic in the four (or fewer)
#' adjacent neighbours in its Rook's neibourhood. This arithmetically averaged
#' difference is then returned under the column header of 'Gradient'. A negative
#' gradient estimate indicates that the statistic in the central cell is higher
#' than that in neighbouring cells whereas a positive gradient estimate
#' indicates the opposite.
#'
#' @param vfdf A data frame as returned by \code{\link{SubgridStats}}
#'
#' @param statistic desired output statistic: It should be one of "mean", "var",
#'   or "sum". Default setting is mean.
#'
#' @return A data frame similar to vfdf except that it includes an additional
#'   column called Gradient as described above.
#'
#' @export
#'
#' @examples
#' # creating pattern patterns
#' Mat1 <- matrix(rep(0,9*9), nrow = 9)
#' Mat1[c(4:6), c(4:6)] <- 2
#' Mat1[c(4:6), c(1:3)] <- 1
#' Mat1[c(1:3), c(4:6)] <- 1
#' Mat1[c(7:9), c(4:6)] <- 1
#' Mat1[c(4:6), c(7:9)] <- 1
#' Mat1
#'
#' Rast1 <- terra::rast(Mat1)
#' terra::plot(Rast1)
#'
#' # calculating the mean in 9 subgrids
#' (statsdf1 <- SubgridStats(Rast1, factv1 = 3, facth1 = 3, statistic = "mean"))
#'
#' # computing the gradient statistic on the mean
#' (graddf1 <- RooksGradient(statsdf1, statistic = "mean"))
#' # the Gradient statistic in the central grid in row 5 should
#' # be equal to negative one
#'
RooksGradient <- function(vfdf, statistic = "mean") {
  if (dim(vfdf)[1] < 5) {
    stop("Number of rows in Vfdf must be at least 5")
  }
  if (is.element(statistic, c("mean", "var", "sum")) == FALSE) {
    stop("statistic must be 'mean', 'var', or 'sum'")
  }

  # creating output data frame
  vfdfout <- vfdf
  vfdfout$Gradient <- rep(0, dim(vfdfout)[1])

  # calculate the spacing between grids
  diffx <- abs(as.numeric(outer(vfdf$colcent, vfdf$colcent, FUN = "-")))
  diffx[diffx == 0] <- NA
  facth <- min(diffx, na.rm = TRUE)
  diffy <- abs(as.numeric(outer(vfdf$rowcent, vfdf$rowcent, FUN = "-")))
  diffy[diffy == 0] <- NA
  factv <- min(diffy, na.rm = TRUE)

  if (factv != facth) {
    stop("Currently, factv must equal facth (see DispField function documentation)")
  }

  # calculate the distance matrix
  distmat <- sqrt((outer(vfdfout$colcent, vfdfout$colcent, FUN = "-")^2) +
                    outer(vfdfout$rowcent, vfdfout$rowcent, FUN = "-") ^2)

  # marking non-neighbours as zero
  diag(distmat) <- 0
  distmat[distmat > factv] <- 0

  # cycling through to compute the gradient statistic
  for (i in 1:dim(vfdfout)[1]) {
    if (statistic == "mean") {
      vfdfout$Gradient[i] = mean(vfdfout$Mean[distmat[i,] > 0] - vfdfout$Mean[i], na.rm = TRUE)
    }

    if (statistic == "var") {
      vfdfout$Gradient[i] = mean(vfdfout$Var[distmat[i,] > 0] - vfdfout$Var[i], na.rm = TRUE)
    }

    if (statistic == "sum") {
      vfdfout$Gradient[i] = mean(vfdfout$Sum[distmat[i,] > 0] - vfdfout$Sum[i], na.rm = TRUE)
    }
  }

  return(vfdfout)
}
