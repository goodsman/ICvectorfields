#' Classify Rook's Neighbours Comprising Spread Patterns in Vector Fields
#'
#' After running the \code{\link{PatternDetect}} function, this function enables
#' classification of neighbour cells. Because the pattern detect function
#' classifies central cells according to the patterns of vector direction in
#' their Rook's neighbourhood, the neighbouring grid locations that comprise the
#' pattern are not labeled. This is remedied by the RooksNeighFind function
#' which classifies neighbour cells around focal grids classified with one of
#' the four patterns that the \code{\link{PatternDetect}} function is able to
#' recognize. The function returns a data frame similar to the input data frame
#' with a column appended. In the appended column, neighbours surrounding focal
#' cells labeled with a particular pattern will be labeled as follows:
#' neighbours of the divergence pattern are labeled with a one, neighbours of
#' the convergence pattern are labeled with a two, neighbours of the partial
#' divergence pattern are labeled with a three, and neighbours of the partial
#' convergence pattern are labeled with a four. In cases where neighbours are
#' shared, the priority order from lowest to highest is four for partial
#' convergence to one for divergence. Thus, a neighbour that is shared between a
#' focal grid classified as a convergence and a nearby focal grid classified as
#' a divergence will be labeled with a one instead of a two.
#'
#' @param vfdf A data frame as returned by \code{\link{PatternDetect}}
#'
#' @return A data frame similar to vfdf except that it includes an additional
#'   column called NeighType as described above.
#'
#' @export
#'
#' @examples
#' # creating convergence/divergence patterns
#' Mat1 <- matrix(rep(0,9*9), nrow = 9)
#' Mat1[3, c(4, 6)] <- 1
#' Mat1[7, c(4, 6)] <- 1
#' Mat1[c(4, 6), 3] <- 1
#' Mat1[c(4, 6), 7] <- 1
#' Mat1
#'
#' Mat2 <- matrix(rep(0,9*9), nrow = 9)
#' Mat2[2, c(4, 6)] <- 1
#' Mat2[8, c(4, 6)] <- 1
#' Mat2[c(4, 6), 2] <- 1
#' Mat2[c(4, 6), 8] <- 1
#' Mat2
#'
#' # rasterizing
#' rast1 <- terra::rast(Mat1)
#' terra::plot(rast1)
#' rast2 <- terra::rast(Mat2)
#' terra::plot(rast2)
#'
#' # Detecting a divergence
#' (VFdf1 <- DispField(rast1, rast2, factv1 = 3, facth1 = 3, restricted = TRUE))
#' (patdf1 <- PatternDetect(VFdf1))
#' (neighdf1 <- RooksNeighFind(patdf1))
#' # Rook's neighbour grids are labeled with a one.
#'
#' # Detecting a convergence
#' (VFdf2 <- DispField(rast2, rast1, factv1 = 3, facth1 = 3, restricted = TRUE))
#' (patdf2 <- PatternDetect(VFdf2))
#' (neighdf2 <- RooksNeighFind(patdf2))
#' # Rook's neighbour grids are labeled with a two.
#'
RooksNeighFind <- function(vfdf) {
  if (is.element("PatternCt", colnames(vfdf)) == FALSE) {
    stop("Must run PatternDetect function first to calculate PatternCt variable")
  }
  if (sum(vfdf$PatternCt, na.rm = TRUE) < 1) {
    stop("Must be at least one candidate (one) in PatternCt column")
  }

  # creating output data frame
  vfdfout <- vfdf
  vfdfout$NeighType <- rep(0, dim(vfdfout)[1])

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

  diag(distmat) <- 0
  distmat[distmat > factv] <- 0

  for (i in 1:dim(vfdfout)[1]) {
    # neighbours of partial convergence pattern
    if (is.na(vfdfout$Pattern[i]) == FALSE & vfdfout$Pattern[i] == "partconv") {
      vfdfout$NeighType[distmat[i,] > 0] <- 4
    }

    # neighbours of partial divergence pattern
    if (is.na(vfdfout$Pattern[i]) == FALSE & vfdfout$Pattern[i] == "partdiv") {
      vfdfout$NeighType[distmat[i,] > 0] <- 3
    }

    # neighbours of convergence pattern
    if (is.na(vfdfout$Pattern[i]) == FALSE & vfdfout$Pattern[i] == "convergence") {
      vfdfout$NeighType[distmat[i,] > 0] <- 2
    }

    # neighbours of divergence pattern
    if (is.na(vfdfout$Pattern[i]) == FALSE & vfdfout$Pattern[i] == "divergence") {
      vfdfout$NeighType[distmat[i,] > 0] <- 1
    }
  }

  return(vfdfout)
}
