# Here are some utility functions, mostly related to matrix manipulation.
# These are called by the main functions in the ICvectorfields package.

## =============================================================================
# This function creates a reflection matrix with ones along the diagonal that
# traverses from the lower left corner to the upper right corner of the square
# matrix of dimension dim1 X dim1.
ReflMat <- function(dim1) {
  Rmat <- matrix(rep(0, dim1^2), nrow = dim1)
  rowmat <- matrix(rep(1:dim1, dim1), nrow = dim1, byrow = F)
  colmat <- matrix(rep(1:dim1, dim1), nrow = dim1, byrow = T)
  Rmat[rowmat + colmat == (dim1 + 1)] <- 1
  return(Rmat)
}

## =============================================================================
# It is much faster to do the matrix flipping using linear algebra
# This flips horizontally and vertically along the central axis.
# Currently the matrix must be square.
FlipMat <- function(inputmat) {
  flipM <- ReflMat(dim1 = dim(inputmat)[1])
  outputmat <- flipM %*% inputmat %*% flipM
  return(outputmat)
}

## =============================================================================
# This function pads any shape of matrix to make a square
# matrix by surrounding the input matrix with a pad of zeros.
# It makes sure that the matrix has an even number of rows and
# columns which facilitates reflections when computing cross-
# covariance.
PadMat <- function(inputmat) {

  # Adding zeros to the left or top to make a square matrix
  if (dim(inputmat)[1] > dim(inputmat)[2]) {
    LeftPad <- matrix(rep(0, (dim(inputmat)[1] - dim(inputmat)[2]) * dim(inputmat)[1]),
      nrow = dim(inputmat)[1]
    )

    inputmat <- cbind(LeftPad, inputmat)
  }
  if (dim(inputmat)[2] > dim(inputmat)[1]) {
    TopPad <- matrix(rep(0, (dim(inputmat)[2] - dim(inputmat)[1]) * dim(inputmat)[2]),
      nrow = (dim(inputmat)[2] - dim(inputmat)[1])
    )

    inputmat <- rbind(TopPad, inputmat)
  }

  # Adding extra zeros
  dimnew <- dim(inputmat)[1]
  PadMat <- matrix(rep(0, (3 * dimnew)^2), nrow = 3 * dimnew)
  PadMat[(dimnew + 1):(2 * dimnew), (dimnew + 1):(2 * dimnew)] <- inputmat

  # unsuring an even number of rows and columns
  if (round(dim(PadMat)[1] / 2.0) != (dim(PadMat)[1] / 2.0)) {
    PadMat <- cbind(PadMat, matrix(rep(0, dim(PadMat)[1]), nrow = dim(PadMat)[1]))
    PadMat <- rbind(PadMat, matrix(rep(0, dim(PadMat)[2]), ncol = dim(PadMat)[2]))
  }

  return(PadMat)
}

## =============================================================================
# This function takes an imput matrix and bounding box coordinates and returns
# a matrix of the same size with zeros at all locations outside the bounding box
ExtractMat <- function(inputmat, rowmin, rowmax, colmin, colmax) {
  inputsub <- inputmat[rowmin:rowmax, colmin:colmax]
  Outmat <- matrix(rep(0, dim(inputmat)[1] * dim(inputmat)[2]),
    nrow = dim(inputmat)[1], ncol = dim(inputmat)[2]
  )
  Outmat[rowmin:rowmax, colmin:colmax] <- inputsub
  return(Outmat)
}

## =============================================================================
# This function takes as input a matrix and returns a data-frame object
# containing references to cells in the original matrix by row and column.
# It does so by subsetting the input matrix into smaller matrices of dimension
# factv X facth, where factv and facth must be odd integers.
# It starts in the upper left corner and proceeds until it hits one of
# the boundary edges. The output data frame contains row and column indices
# referring to the edges of each sub-matrix of dimension factv X facth as well
# as the central cell in the following order: central cell row and column, focal
# cell row minimum and row maximum, focal cell column minimum and column maximum
ThinMat <- function(inputmat, factv, facth) {
  if (factv / 2 == round(factv / 2)) {
    stop("factv and facth must be odd integers")
  }
  if (facth / 2 == round(facth / 2)) {
    stop("factv and facth must be odd integers")
  }
  if (round(factv) != factv || round(facth) != facth) {
    stop("factv and facth must be integers")
  }

  # First finding all centroid references
  dx <- floor(facth / 2) + 1
  dy <- floor(factv / 2) + 1
  Outdf <- expand.grid(
    seq(dy, dim(inputmat)[1], by = (2 * floor(facth / 2) + 1)),
    seq(dx, dim(inputmat)[2], by = (2 * floor(facth / 2) + 1))
  )

  # focal matrix min and max
  colnames(Outdf) <- c("rowcent", "colcent")
  Outdf$frowmin <- Outdf$rowcent - (facth - dy)
  Outdf$frowmax <- Outdf$rowcent + (facth - dy)
  Outdf$fcolmin <- Outdf$colcent - (facth - dx)
  Outdf$fcolmax <- Outdf$colcent + (facth - dx)

  # Removing locations that bump against the edge
  Outdf[Outdf < 1] <- NA
  Outdf <- stats::na.omit(Outdf)
  Outdf <- subset(Outdf, Outdf$frowmax <= dim(inputmat)[1])
  Outdf <- subset(Outdf, Outdf$fcolmax <= dim(inputmat)[2])

  return(Outdf)
}

## =============================================================================
# This function takes a raster file in terra format and converts it to a
# rectangular matrix of the same dimensions as the input raster. It also converts
# NA values to zero.
RastToMatrix <- function(inrast) {
  outmat <- matrix(as.vector(inrast), nrow = dim(inrast)[1], ncol = dim(inrast)[2], byrow = T)
  outmat[is.na(outmat) == TRUE] <- 0
  return(outmat)
}
