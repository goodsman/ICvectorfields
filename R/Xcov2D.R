#' Cross-covariance in two spatial dimensions
#'
#' This function efficiently computes two dimensional cross-covariance of two
#' equal dimensioned matrices of real numbers using efficient discrete fast
#' Fourier trasforms.
#'
#' The algorithm first pads each matrix with zeros so that the outer edges of
#' the matrices do not interact with one another due to the circular nature of
#' the discrete fast Fourier transform. Cross-covariance calculations require
#' computation of the complex conjugate of one of the two imput matrices.
#' Assuming all of it's elements are real, computing the complex conjugate is
#' equivalent to flipping the matrix in the horizontal and vertical directions.
#' Then to compute cross-covariance, the first matrix is convolved with the
#' flipped second matrix as described in the convolution theorem.
#'
#' This function is called by the main functions that compute displacement
#' fields and vector fields and is included here primarily for demonstration
#' purposes. Specifically, the method for computing the magnitude and direction
#' of shifts is demonstrated in the examples.
#'
#' The shift that produces the maximum cross-covariance between the two input
#' matrices can be obtained by finding the row and column indices associated
#' with the maximum cross-covariance. The shift in each direction is obtained by
#' subracting one plus the half the dimension of the output matrix (the same for
#' rows and columns) from the row and column values that are associated with the
#' maximum cross-covariance as demonstrated in the examples below. Note that
#' shifts to the right and up are denoted with positive numbers and shifts to
#' the left and down are denoted by negative numbers. This is contrary to some
#' conventions but efficient for producing vector fields. For more details on
#' cross-covariance see
#' \href{https://en.wikipedia.org/wiki/Cross-correlation}{cross-correlation}.
#'
#' @param mat1 a real valued matrix
#' @param mat2 a real valued matrix of equal dimension to mat1
#'
#' @return a real valued matrix showing cross-covariance in each direction
#' @export
#'
#' @examples
#' matrix(c(1:6, rep(0, 3)), nrow = 3); matrix(c(rep(0, 3), 1:6), nrow = 3)
#' dim(Xcov2D(matrix(c(1:6, rep(0, 3)), nrow = 3),
#'     matrix(c(rep(0, 3), 1:6), nrow = 3)))
#' ICvectorfields::GetRowCol(
#'     which.max(Xcov2D(matrix(c(1:6, rep(0, 3)), nrow = 3),
#'         matrix(c(rep(0, 3), 1:6), nrow = 3))),
#'     dim1 = dim(Xcov2D(matrix(c(1:6, rep(0, 3)), nrow = 3),
#'         matrix(c(rep(0, 3), 1:6), nrow = 3)))[1],
#'     dim2 = dim(Xcov2D(matrix(c(1:6, rep(0, 3)), nrow = 3),
#'         matrix(c(rep(0, 3), 1:6), nrow = 3)))[2]
#'         )
#' # This implies that the shift is 6 - (10/2 + 1) in the vertical
#' # direction and 7 - (10/2 + 1) in the horizonatal direction.
Xcov2D <- function(mat1, mat2) {
  if (dim(mat1)[1] != dim(mat2)[1] || dim(mat1)[2] != dim(mat2)[2]) {
    stop("unequal dimensions")
  }

  # padding the matrices, and making sure
  # they are square with an even number of rows and columns
  mat1pad <- PadMat(mat1)
  mat2pad <- PadMat(mat2)

  # flip one matrix in the x and y directions
  mat2new <- FlipMat(mat2pad)

  # Now I convolve the two
  fmat3 <- fftwtools::fftw_r2c_2d(mat1pad) * fftwtools::fftw_r2c_2d(mat2new)
  mat3 <- fftwtools::fftw2d(fmat3, inverse = 1) / length(mat1pad)^2

  # flipping
  Index1 <- dim(mat1pad)[1] / 2
  LR <- mat3[(Index1 + 1):(2 * Index1), (Index1 + 1):(2 * Index1)]
  LR <- FlipMat(LR)
  UL <- mat3[1:Index1, 1:Index1]
  UL <- FlipMat(UL)
  LL <- mat3[(Index1 + 1):(2 * Index1), 1:Index1]
  LL <- FlipMat(LL)
  UR <- mat3[1:Index1, (Index1 + 1):(2 * Index1)]
  UR <- FlipMat(UR)

  # reassembling
  OutMatlower <- cbind(LL, LR)
  OutMatupper <- cbind(UL, UR)
  OutMat <- rbind(OutMatupper, OutMatlower)

  return(Re(OutMat))
}
