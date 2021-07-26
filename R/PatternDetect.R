#' Detect Patterns in Vector Fields
#'
#' Detect patterns in vector fields represented on a grid by looking in the
#' rook's neighbourhood of each grid cell. Three patterns are detected: convergences
#' occur when the vectors in the four adjacent cells point towards the focal
#' cell; divergences occur when the vectors in the four adjacent cells point
#' away from the focal cell; waves occur when directions of vectors are
#' consistent across all of the adjacent cells as well as the focal cell. In the
#' wave pattern vectors can point up and right, up and left, down and right, or
#' down and left.
#'
#' @param vfdf A data frame as returned by \code{\link{DispField}},
#'   \code{\link{DispFieldST}}, or \code{\link{DispFieldSTall}} with at least
#'   five rows (more is better)
#'
#' @return A data frame as returned by \code{\link{DispField}},
#'   \code{\link{DispFieldST}}, or \code{\link{DispFielSTall}}, with an
#'   additional column called Pattern in which the patterns around each focal
#'   cell are categorized as convergence, divergence, wave or NA
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
#'
#' # Detecting a convergence
#' (VFdf2 <- DispField(rast2, rast1, factv1 = 3, facth1 = 3, restricted = TRUE))
#' (patdf2 <- PatternDetect(VFdf2))
PatternDetect <- function(vfdf) {
  if(dim(vfdf)[1] < 5) {
    stop("vfdf data frame must have at least 5 rows")
  }

  # generating a output column for pattern analysis
  vfdfout <- vfdf
  vfdfout$Pattern <- rep(NA, dim(vfdf)[1])

  # computing the distance between centres of grids
  diffx <- diff(vfdf$colcent)
  diffx[diffx <= 0] <- NA
  facth <- min(diffx, na.rm = TRUE)
  diffy <- diff(vfdf$rowcent)
  diffy[diffy <= 0] <- NA
  factv <- min(diffy, na.rm = TRUE)

  for (i in 1:dim(vfdf)[1]) {
    # Computing distance between focal cell and all other cells
    DistVec <- sqrt(((vfdf$rowcent - vfdf$rowcent[i])^2) + (vfdf$colcent - vfdf$colcent[i])^2)

    # finding displacement in x and y in the cell above
    Updx <- vfdf$dispx[is.na(match(vfdf$rowcent, vfdf$rowcent[i] - factv)) == FALSE &
                         DistVec == factv]
    Updy <- vfdf$dispy[is.na(match(vfdf$rowcent, vfdf$rowcent[i] - factv)) == FALSE &
                         DistVec == factv]

    # finding displacement in x and y in the cell below
    Downdx <- vfdf$dispx[is.na(match(vfdf$rowcent, vfdf$rowcent[i] + factv)) == FALSE &
                           DistVec == factv]
    Downdy <- vfdf$dispy[is.na(match(vfdf$rowcent, vfdf$rowcent[i] + factv)) == FALSE &
                           DistVec == factv]

    # finding displacement in x and y in the cell to the left
    Leftdx <- vfdf$dispx[is.na(match(vfdf$colcent, vfdf$colcent[i] - facth)) == FALSE &
                           DistVec == factv]
    Leftdy <- vfdf$dispy[is.na(match(vfdf$colcent, vfdf$colcent[i] - facth)) == FALSE &
                           DistVec == factv]

    # finding displacement in x and y in the cell to the right
    Rightdx <- vfdf$dispx[is.na(match(vfdf$colcent, vfdf$colcent[i] + facth)) == FALSE &
                            DistVec == factv]
    Rightdy <- vfdf$dispy[is.na(match(vfdf$colcent, vfdf$colcent[i] + facth)) == FALSE &
                            DistVec == factv]

    if (length(Updx) == 1 & length(Updy) == 1 &
        length(Downdx) == 1 & length(Downdy) == 1 &
        length(Leftdx) == 1 & length(Leftdy) == 1 &
        length(Rightdx) == 1 & length(Rightdy) == 1) {
      # Finding convergences
      if (Updy < 0 & Downdy > 0 & Leftdx > 0 & Rightdx < 0) {
        vfdfout$Pattern[i] <- "convergence"
      }

      # Finding divergences
      if (Updy > 0 & Downdy < 0 & Leftdx < 0 & Rightdx > 0) {
        vfdfout$Pattern[i] <- "divergence"
      }

      # Up-right waves
      if (Updy > 0 & Downdy > 0 & Leftdy > 0 & Rightdy > 0 & vfdf$dispy[i] > 0 &
          Updx > 0 & Downdx > 0 & Leftdx > 0 & Rightdx > 0 & vfdf$dispx[i] > 0) {
        vfdfout$Pattern[i] <- "wave"
      }

      # Up-left waves
      if (Updy > 0 & Downdy > 0 & Leftdy > 0 & Rightdy > 0 & vfdf$dispy[i] > 0 &
          Updx < 0 & Downdx < 0 & Leftdx < 0 & Rightdx < 0 & vfdf$dispx[i] < 0) {
        vfdfout$Pattern[i] <- "wave"
      }

      # Down-right waves
      if (Updy < 0 & Downdy < 0 & Leftdy < 0 & Rightdy < 0 & vfdf$dispy[i] < 0 &
          Updx > 0 & Downdx > 0 & Leftdx > 0 & Rightdx > 0 & vfdf$dispx[i] > 0) {
        vfdfout$Pattern[i] <- "wave"
      }

      # Down-left waves
      if (Updy < 0 & Downdy < 0 & Leftdy < 0 & Rightdy < 0 & vfdf$dispy[i] < 0 &
          Updx < 0 & Downdx < 0 & Leftdx < 0 & Rightdx < 0 & vfdf$dispx[i] < 0) {
        vfdfout$Pattern[i] <- "wave"
      }
    }
  }

  return(vfdfout)
}
