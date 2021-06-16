#' Creating a raster stack from formatted datasets
#'
#' This function converts the data that accompany the ICvectorfields R package
#' to a raster stack. The raster stack is the only accepted data input format
#' for the following ICvectorfields functions: \code{\link{DispFieldST}},
#' \code{\link{DispFieldSTbb}}, \code{\link{DispFieldSTall}},
#' \code{\link{DispFieldSTbball}}.
#'
#' Once a raster stack has been created, individual layers can be subsetted
#' using rasterstack[\[index\]], where index is an integer index for the third
#' dimension of the raster stack.
#'
#' @param inputdf a data frame object in which the first column is longitude (or
#'   x coordinate), the second column is latitude (or y coordinate), and all of
#'   the subsequent columns represent a measure of population abundance or
#'   density at a unique instance of time. Each row of the input data frame,
#'   therefore, represents a unique spatial location, which should be on an
#'   evenly spaced grid. Note, however, that not all grid locations need to have
#'   observations; some grid locations can have values of NA or can be missing
#'   entirely.
#'
#' @return The function returns a raster stack constructed using inputdf. Each
#'   layer in the stack corresponds to a column of the input dataset (after the
#'   first two columns, which are longitude and latitude). The extent of all of
#'   the rasters in the stack is constructed using the minimum and maximum
#'   longitudes and latitudes.
#' @export
#'
#' @examples
#'
#' # creating random data in the correct data format
#' xyzdf <- expand.grid(x = c(1:3), y = c(1:3))
#' xyzdf$z1 <- runif(9)
#' xyzdf$z2 <- runif(9)
#' xyzdf$z3 <- runif(9)
#'
#' zstack <- RastStackData(xyzdf)
#'
#' dim(zstack)
#' terra::plot(zstack[[1]])
#' terra::plot(zstack[[2]])
#' terra::plot(zstack[[3]])
RastStackData <- function(inputdf) {

  # initializing
  inputdfstack <- c()

  # creating the stack
  for (i in 3:dim(inputdf)[2]) {
    inputdfrast <- terra::rast(as.matrix(inputdf[, c(1:2, i)]), type = "xyz")
    terra::ext(inputdfrast) <- c(min(inputdf[, 1]), max(inputdf[, 1]), min(inputdf[, 2]), max(inputdf[, 2]))
    inputdfstack <- c(inputdfstack, inputdfrast)
  }

  # ensuring the terra recognizes output as rast
  inputdfstack <- terra::rast(inputdfstack)

  return(inputdfstack)
}
