#' Larch budmoth defoliation data from the ncf R package
#'
#' These are the data in Bjornstad et al. (2002). The data
#' are formatted such that they can easily be converted to
#' a raster stack using
#' ICvectorfields::RastStackData(LBMfromncfpkg).
#'
#' @usage data(LBMfromncfpkg)
#' @usage ICvectorfields::RastStackData(LBMfromncfpkg)
#'
#' @format A data-frame with 135 rows and 40 columns. The first two are the x-
#'   and y-coordinates (in m), the following 38 columns represent defoliation in
#'   years 1961 through 1998.
#'
#' @source \url{https://cran.r-project.org/web/packages/ncf/index.html}
#'
#' @references Bjornstad, O.N., Peltonen, M., Liebhold, A.M., and Baltensweiler,
#'   W. (2002) Waves of larch budmoth outbreaks in the European Alps. Science,
#'   298, 1020-1023. \url{https://doi.org/10.1126/science.1075182}
#'
#' @keywords datasets

"LBMfromncfpkg"

#' Copepod abundance data
#'
#' Data are an amalgamation of many years of sampling. They show seasonal
#' patterns of Copepod abundance in the world's oceans. The original dataset was
#' clipped so that the region around Japan is central in the dataset supplied
#' with the ICvectorfields package. Sampling grid is 0.25 X 0.25 degrees. The
#' data are formatted such that they can easily be converted to a raster stack
#' using ICvectorfields::RastStackData(Copepod).
#'
#' @usage data(Copepod)
#' @usage ICvectorfields::RastStackData(Copepod)
#'
#' @format A data-frame with 321 rows and 5 columns.
#' \describe{
#'   \item{Longitude}{in degrees}
#'   \item{Latitude}{in degrees}
#'   \item{wtmass5}{wet mass in May in milligrams per cubic meter}
#'   \item{wtmass6}{wet mass in June in milligrams per cubic meter}
#'   \item{wtmass7}{wet mass in July in milligrams per cubic meter}
#' }
#'
#' @source \url{https://www.st.nmfs.noaa.gov/copepod/about/spatial-fields.html}
#'
#' @references Moriarty, R., and O'Brien,  T. D. (2013) Distribution of
#' mesozooplankton biomass in the global ocean. Earth Syst. Sci. Data,
#' 5, 45-55. \url{https://doi.org/10.5194/essd-5-45-2013}
#'
#' @keywords datasets

"Copepod"

#' Simulated movement data
#'
#' Data based on a partial differential equation were simulated using the
#' ReacTran R package. The simulation algorith uses a finite differencing scheme
#' with backwards differencing. The model used for simulation is a reaction
#' diffusion-advection equation in which the advection term is variable in space
#' but diffusion and reactions are constant in space:
#'
#' \deqn{\frac{\partial(c)}{\partial{t}} = \nabla (D \nabla c) - \nabla (v c) + r}
#'
#' The parameters used in the partial differntial equation above are
#'
#' D = 0.01 per squared spatial unit
#'
#' r = 0.5 per unit time
#'
#' v (advection is variable in space): in the upper left quadrant of the
#' square domain v = (0.2, 0); in the upper right quadrant v = (0, -0.2);
#' in the lower right quadrant v = (-0.2, 0); in the lower right quadrant
#' v = (0, 0.2). Obviously v is discontinous at the quadrant boundaries,
#' which causes some interesting model behaviour that is limited by
#' considering only the first six time steps such that the bulk of the
#' concentration in each quadrant does not cross a quadrant boundary.
#'
#' The intial condition at time = 0 is a concentration of one unit per
#' arbitrary unit of volume in the central cell of each quadrant.
#'
#' External boundary conditions are zero-gradient (reflecting).
#'
#' The data are formatted such that they can easily be converted to a raster
#' stack using ICvectorfields::RastStackData(SimData).
#'
#' @usage data(SimData)
#' @usage ICvectorfields::RastStackData(SimData)
#'
#' @format A data-frame with 40804 rows and 8 columns.
#' \describe{
#'   \item{xcoord}{in arbitrary units}
#'   \item{ycoord}{in arbitrary units}
#'   \item{t1}{concentration in arbitrary units at t = 1}
#'   \item{t2}{concentration in arbitrary units at t = 2}
#'   \item{t3}{concentration in arbitrary units at t = 3}
#'   \item{t4}{concentration in arbitrary units at t = 4}
#'   \item{t5}{concentration in arbitrary units at t = 5}
#'   \item{t6}{concentration in arbitrary units at t = 6}
#' }
#'
#' @keywords datasets

"SimData"
