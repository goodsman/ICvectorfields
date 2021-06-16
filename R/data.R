#' Simulated movement data
#'
#' Data based on a partial differential equation were simulated using the
#' ReacTran R package (see details).
#'
#' The simulation algorithm uses a finite differencing scheme with backwards
#' differencing. The model used for simulation is a reaction diffusion-advection
#' equation in which the advection term is variable in space but diffusion and
#' reactions are constant in space see
#' \href{https://en.wikipedia.org/wiki/Convection%E2%80%93diffusion_equation}{convection-diffusion
#' equation} for an example.
#'
#' The parameters used in the general partial differntial equation in the link
#' above are
#'
#' D = 0.01 per squared spatial unit
#'
#' R = 0.5 per unit time
#'
#' v (advection is variable in space): in the upper left quadrant of the square
#' domain v = (0.2, 0); in the upper right quadrant v = (0, -0.2); in the lower
#' right quadrant v = (-0.2, 0); in the lower right quadrant v = (0, 0.2).
#' Obviously v is discontinous at the quadrant boundaries, which causes some
#' interesting model behaviour that is limited by considering only the first six
#' time steps such that the bulk of the concentration in each quadrant does not
#' cross a quadrant boundary.
#'
#' The intial condition at time = 0 is a concentration of one unit per arbitrary
#' unit of volume in the central cell of each quadrant.
#'
#' External boundary conditions are zero-gradient (reflecting).
#'
#' The data are formatted such that they can easily be converted to a raster
#' stack using ICvectorfields::RastStackData(SimData).
#'
#' @usage data(SimData)
#'
#' @format A data-frame with 40804 rows and 8 columns. \describe{
#'   \item{xcoord}{in arbitrary units} \item{ycoord}{in arbitrary units}
#'   \item{t1}{concentration in arbitrary units at t = 1}
#'   \item{t2}{concentration in arbitrary units at t = 2}
#'   \item{t3}{concentration in arbitrary units at t = 3}
#'   \item{t4}{concentration in arbitrary units at t = 4}
#'   \item{t5}{concentration in arbitrary units at t = 5}
#'   \item{t6}{concentration in arbitrary units at t = 6} }
#'
#' @keywords datasets

"SimData"
