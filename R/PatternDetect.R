#' Detect Patterns in Vector Fields
#'
#' Detect patterns in vector fields represented on a grid by looking in the
#' rook's neighbourhood of each grid cell. Four patterns are detected:
#' convergences occur when the vectors in the four adjacent cells in the rook's
#' neighbourhood point towards the focal cell; divergences occur when the
#' vectors in the four adjacent cells point away from the focal cell; Partial
#' convergences occur when three of the four vectors point towards the focal
#' cell and the final vector points neither towards nor away from the focal
#' cell; Partial divergences occur when three of the four vectors point away the
#' focal grid cell and the final vector points neither towards nor away from the
#' focal grid. For all of the patterns above a sub-pattern is specified if all
#' arrows point clockwise or counter-clockwise.
#'
#' @param vfdf A data frame as returned by \code{\link{DispField}},
#'   \code{\link{DispFieldST}}, or \code{\link{DispFieldSTall}} with at least
#'   five rows (more is better)
#'
#' @return A data frame as returned by \code{\link{DispField}},
#'   \code{\link{DispFieldST}}, or \code{\link{DispFieldSTall}}, with three
#'   additional columns. The first additional column is called Pattern in which
#'   the patterns around each focal cell are categorized as convergence,
#'   divergence, partial convergence, partial divergence, or NA. The second
#'   additional column, called SubPattern, indicates whether all arrows point
#'   clockwise or counter-clockwise. The third additional column is called
#'   PatternCt, which contains a one if all four neighbourhood grid cells
#'   contain displacement estimates, and a NA otherwise.
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

      # Finding convergences
      if (Updy < 0 & Downdy > 0 & Leftdx > 0 & Rightdx < 0) {
        vfdfout$Pattern[i] <- "convergence"
      }

      # Finding partial convergence type 1
      if (Updy == 0 & Downdy > 0 & Leftdx > 0 & Rightdx < 0) {
        vfdfout$Pattern[i] <- "partconv"
      }

      # Finding partial convergences type 2
      if (Updy < 0 & Downdy == 0 & Leftdx > 0 & Rightdx < 0) {
        vfdfout$Pattern[i] <- "partconv"
      }

      # Finding partial convergences type 3
      if (Updy < 0 & Downdy > 0 & Leftdx == 0 & Rightdx < 0) {
        vfdfout$Pattern[i] <- "partconv"
      }

      # Finding partial convergences type 4
      if (Updy < 0 & Downdy > 0 & Leftdx > 0 & Rightdx == 0) {
        vfdfout$Pattern[i] <- "partconv"
      }

      # Finding divergences
      if (Updy > 0 & Downdy < 0 & Leftdx < 0 & Rightdx > 0) {
        vfdfout$Pattern[i] <- "divergence"
      }

      # Finding partial divergences type 1
      if (Updy == 0 & Downdy < 0 & Leftdx < 0 & Rightdx > 0) {
        vfdfout$Pattern[i] <- "partdiv"
      }

      # Finding partial divergences type 2
      if (Updy > 0 & Downdy == 0 & Leftdx < 0 & Rightdx > 0) {
        vfdfout$Pattern[i] <- "partdiv"
      }

      # Finding partial divergences type 3
      if (Updy > 0 & Downdy < 0 & Leftdx == 0 & Rightdx > 0) {
        vfdfout$Pattern[i] <- "partdiv"
      }

      # Finding partial divergences type 4
      if (Updy > 0 & Downdy < 0 & Leftdx < 0 & Rightdx == 0) {
        vfdfout$Pattern[i] <- "partdiv"
      }

      # Finding sub-patterns that are clockwise in all four neighbours
      if (Updx > 0 & Downdx < 0 & Leftdy > 0 & Rightdy < 0) {
        vfdfout$SubPattern[i] <- "clockwise"
      }

      # Finding sub-patterns that are counter-clockwise in all four neighbours
      if (Updx < 0 & Downdx > 0 & Leftdy < 0 & Rightdy > 0) {
        vfdfout$SubPattern[i] <- "counter"
      }

    }
  }

  return(vfdfout)
}
