#' Xcov2D
#'
#' This function efficiently computes 2D cross-covariance of two equal
#' dimensioned matrices of real numbers using efficient discrete fast Fourier
#' trasforms.
#'
#' The algorithm first pads each matrix so that the outer edges of the matrices
#' do not interact with one another due to the circular nature of the discrete
#' fast Fourier transform. Cross-correlation is then computed by computing the
#' complex conjugate of the second matrix assuming all of it's elements are
#' real. This operation is equivalent to flipping the matrix in the horizontal
#' and vertical directions. Then the first matrix is convolved with the flipped
#' second matrix using the convolution theorem.
#'
#' @param mat1 a numeric real valued matrix
#' @param mat2 a numeric real valued matrix of equal dimension to mat1
#'
#' @return a real valued matrix showing cross-covariance in each direction
#' @export
#'
#' @examples
#' dim(Xcov2D(matrix(c(1:6, rep(0, 3)), nrow = 3), matrix(c(rep(0, 3), 1:6), nrow = 3)))
#' GetRowCol(which.max(Xcov2D(matrix(c(1:6, rep(0, 3)), nrow = 3), matrix(c(rep(0, 3), 1:6), nrow = 3))),
#'     dim1 = dim(Xcov2D(matrix(c(1:6, rep(0, 3)), nrow = 3), matrix(c(rep(0, 3), 1:6), nrow = 3)))[1],
#'     dim2 = dim(Xcov2D(matrix(c(1:6, rep(0, 3)), nrow = 3), matrix(c(rep(0, 3), 1:6), nrow = 3)))[2])
#' max(Xcov2D(matrix(c(1:6, rep(0, 3)), nrow = 3), matrix(c(rep(0, 3), 1:6), nrow = 3)))
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
