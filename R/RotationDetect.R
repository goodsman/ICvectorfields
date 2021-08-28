#' Detect Rotating Patterns in Vector Fields
#'
#' Detect patterns in vector fields represented on a grid by looking in the
#' rook's neighbourhood of each grid cell. This function is analogous to
#' \code{\link{PatternDetect}}, except that it detects rotational patterns. Four
#' patterns are detected: clockwise rotation when rotation in all four
#' neighbbour grids appears clockwise, counter clockwise rotation when rotation
#' in all four neighbour grids appears counter-clockwise, and partial clockwise
#' and counter-clockwise rotation, when all but one of the four adjacent
#' neighbour cells has vectors that indicate rotation. For all of the patterns
#' above a sub-pattern is specified as convergence when all of the vectors in
#' the four adjacent grids point towards the focal cell or a divergence when all
#' of the vectors in the four adjacent grids point away from the focal cell.
#'
#' @param vfdf A data frame as returned by \code{\link{DispField}},
#'   \code{\link{DispFieldST}}, or \code{\link{DispFieldSTall}} with at least
#'   five rows (more is better)
#'
#' @return A data frame as returned by \code{\link{DispField}},
#'   \code{\link{DispFieldST}}, or \code{\link{DispFieldSTall}}, with three
#'   additional columns. The first additional column is called Pattern in which
#'   the patterns around each focal cell are categorized as clockwise,
#'   counter-clockwise, partial clockwise, partial counter-clockwise, or NA. The
#'   second additional column, called SubPattern, indicates whether all arrows
#'   point towards (convergence) or away (divergence) from the focal cell. The
#'   third additional column is called PatternCt, which contains a one if all
#'   four neighbourhood grid cells contain displacement estimates, and a NA
#'   otherwise.
#' @export
#'
#' @examples
#' # creating rotation patterns
#' Mat1 <- matrix(rep(0,9*9), nrow = 9)
#' Mat1[c(1:3), 4] <- 1
#' Mat1[c(7:9), 6] <- 1
#' Mat1[4, c(7:9)] <- 1
#' Mat1[6, c(1:3)] <- 1
#' Mat1
#'
#' Mat2 <- matrix(rep(0,9*9), nrow = 9)
#' Mat2[c(1:3), 5] <- 1
#' Mat2[c(7:9), 5] <- 1
#' Mat2[5, c(7:9)] <- 1
#' Mat2[5, c(1:3)] <- 1
#' Mat2
#'
#' # rasterizing
#' rast1 <- terra::rast(Mat1)
#' terra::plot(rast1)
#' rast2 <- terra::rast(Mat2)
#' terra::plot(rast2)
#'
#' # Detecting clockwise rotation
#' (VFdf1 <- DispField(rast1, rast2, factv1 = 3, facth1 = 3, restricted = TRUE))
#' (patdf1 <- RotationDetect(VFdf1))
#'
#' # Detecting counter-clockwise rotation
#' (VFdf2 <- DispField(rast2, rast1, factv1 = 3, facth1 = 3, restricted = TRUE))
#' (patdf2 <- RotationDetect(VFdf2))
RotationDetect <- function(vfdf) {
  if(dim(vfdf)[1] < 5) {
    stop("vfdf data frame must have at least 5 rows")
  }

  # generating a output column for pattern analysis
  vfdfout <- vfdf
  vfdfout$Pattern <- rep(NA, dim(vfdf)[1])
  vfdfout$SubPattern <- rep(NA, dim(vfdf)[1])
  vfdfout$PatternCt <- rep(NA, dim(vfdf)[1])
  vfdfout$dispx[is.na(vfdfout$dispx) == TRUE] = 0.0
  vfdfout$dispy[is.na(vfdfout$dispy) == TRUE] = 0.0

  # computing the distance between centres of grids
  diffx <- abs(as.numeric(outer(vfdf$colcent, vfdf$colcent, FUN = "-")))
  diffx[diffx == 0] <- NA
  facth <- min(diffx, na.rm = TRUE)
  diffy <- abs(as.numeric(outer(vfdf$rowcent, vfdf$rowcent, FUN = "-")))
  diffy[diffy == 0] <- NA
  factv <- min(diffy, na.rm = TRUE)

  for (i in 1:dim(vfdfout)[1]) {
    # Computing distance between focal cell and all other cells
    DistVec <- sqrt(((vfdfout$rowcent - vfdfout$rowcent[i])^2) + (vfdfout$colcent - vfdfout$colcent[i])^2)

    # finding displacement in x and y in the cell above
    Updx <- vfdfout$dispx[is.na(match(vfdfout$rowcent, vfdfout$rowcent[i] - factv)) == FALSE &
                            DistVec == factv]
    Updy <- vfdfout$dispy[is.na(match(vfdfout$rowcent, vfdfout$rowcent[i] - factv)) == FALSE &
                            DistVec == factv]

    # finding displacement in x and y in the cell below
    Downdx <- vfdfout$dispx[is.na(match(vfdfout$rowcent, vfdfout$rowcent[i] + factv)) == FALSE &
                              DistVec == factv]
    Downdy <- vfdfout$dispy[is.na(match(vfdfout$rowcent, vfdfout$rowcent[i] + factv)) == FALSE &
                              DistVec == factv]

    # finding displacement in x and y in the cell to the left
    Leftdx <- vfdfout$dispx[is.na(match(vfdfout$colcent, vfdfout$colcent[i] - facth)) == FALSE &
                              DistVec == factv]
    Leftdy <- vfdfout$dispy[is.na(match(vfdfout$colcent, vfdfout$colcent[i] - facth)) == FALSE &
                              DistVec == factv]

    # finding displacement in x and y in the cell to the right
    Rightdx <- vfdfout$dispx[is.na(match(vfdfout$colcent, vfdfout$colcent[i] + facth)) == FALSE &
                               DistVec == factv]
    Rightdy <- vfdfout$dispy[is.na(match(vfdfout$colcent, vfdfout$colcent[i] + facth)) == FALSE &
                               DistVec == factv]

    # calculating speed in each of the neighboring cells
    if (length(Updx) == 1 & length(Updy) == 1) {
      UpSpeed = sqrt((Updx^2) + Updy^2)
    } else {
      UpSpeed = 0.0
    }
    if (length(Downdx) == 1 & length(Downdy) == 1) {
      DownSpeed = sqrt((Downdx^2) + Downdy^2)
    } else {
      DownSpeed = 0.0
    }
    if (length(Leftdx) == 1 & length(Leftdy) == 1) {
      LeftSpeed = sqrt((Leftdx^2) + Leftdy^2)
    } else {
      LeftSpeed = 0.0
    }
    if (length(Rightdx) == 1 & length(Rightdy) == 1) {
      RightSpeed = sqrt((Rightdx^2) + Rightdy^2)
    } else {
      RightSpeed = 0.0
    }

    if (UpSpeed != 0.0 & DownSpeed != 0.0 &
        LeftSpeed != 0.0 & RightSpeed != 0.0) {

      # indicating whether all four neighbourhood cells have displacement estimates
      vfdfout$PatternCt[i] = 1

      # Finding sub-patterns that are clockwise in all four neighbours
      if (Updx > 0 & Downdx < 0 & Leftdy > 0 & Rightdy < 0) {
        vfdfout$Pattern[i] <- "clockwise"
      }

      # Finding sub-patterns that are clockwise in 3 neighbours type 1
      if (Updx <= 0 & Downdx < 0 & Leftdy > 0 & Rightdy < 0) {
        vfdfout$Pattern[i] <- "partclock"
      }

      # Finding sub-patterns that are clockwise in 3 neighbours type 2
      if (Updx > 0 & Downdx >= 0 & Leftdy > 0 & Rightdy < 0) {
        vfdfout$Pattern[i] <- "partclock"
      }

      # Finding sub-patterns that are clockwise in 3 neighbours type 3
      if (Updx > 0 & Downdx < 0 & Leftdy <= 0 & Rightdy < 0) {
        vfdfout$Pattern[i] <- "partclock"
      }

      # Finding sub-patterns that are clockwise in 3 neighbours type 4
      if (Updx > 0 & Downdx < 0 & Leftdy > 0 & Rightdy >= 0) {
        vfdfout$Pattern[i] <- "partclock"
      }

      # Finding sub-patterns that are counter-clockwise in all four neighbours
      if (Updx < 0 & Downdx > 0 & Leftdy < 0 & Rightdy > 0) {
        vfdfout$Pattern[i] <- "counter"
      }

      # Finding sub-patterns that are counter-clockwise in 3 neighbours type 1
      if (Updx >= 0 & Downdx > 0 & Leftdy < 0 & Rightdy > 0) {
        vfdfout$Pattern[i] <- "partcounter"
      }

      # Finding sub-patterns that are counter-clockwise in 3 neighbours type 2
      if (Updx < 0 & Downdx <= 0 & Leftdy < 0 & Rightdy > 0) {
        vfdfout$Pattern[i] <- "partcounter"
      }

      # Finding sub-patterns that are counter-clockwise in 3 neighbours type 3
      if (Updx < 0 & Downdx > 0 & Leftdy >= 0 & Rightdy > 0) {
        vfdfout$Pattern[i] <- "partcounter"
      }

      # Finding sub-patterns that are counter-clockwise in 3 neighbours type 4
      if (Updx < 0 & Downdx > 0 & Leftdy < 0 & Rightdy <= 0) {
        vfdfout$Pattern[i] <- "partcounter"
      }

      # Finding convergences
      if (Updy < 0 & Downdy > 0 & Leftdx > 0 & Rightdx < 0) {
        vfdfout$SubPattern[i] <- "convergence"
      }

      # Finding divergences
      if (Updy > 0 & Downdy < 0 & Leftdx < 0 & Rightdx > 0) {
        vfdfout$SubPattern[i] <- "divergence"
      }

    }
  }

  return(vfdfout)
}
