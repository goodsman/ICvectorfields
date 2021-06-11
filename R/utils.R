# Here are some utility functions, mostly related to matrix manipulation.
# These are called by the main functions in the ICvectorfields package.

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
