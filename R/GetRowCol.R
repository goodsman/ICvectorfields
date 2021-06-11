#' GetRowCol
#'
#' Here is a function that will find the row and column associated with a vector index.
#'
#' Often when applying functions like the R function which.max(matrix) to a matrix a vector index
#' is returned when the coder would prefer to have a row and index number. This function converts
#' the index to row and column numbers.
#'
#' The function assumes that the elements of the matrix are filled by row (byrow = TRUE),
#' which is the default
#'
#' @param Index an integer vector index
#' @param dim1 integer row dimension of the matrix from which the row and column indices are to be extracted
#' @param dim2 integer column dimension of the matrix from which the row and column indices are to be extracted
#'
#' @return a numeric vector of length two with two integers indicating row and column respectively
#' @export
#'
#' @examples
#' GetRowCol(6, dim1 = 3, dim2 = 3) # should return c(3, 2)
GetRowCol <- function(Index, dim1, dim2) {
  Col1 <- ceiling(Index / dim1)
  Row1 <- Index - (Col1 - 1) * dim1
  return(c(Row1, Col1))
}
