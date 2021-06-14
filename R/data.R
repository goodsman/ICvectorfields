#' Larch budmoth defoliation data from the ncf R package
#'
#' These are the data in Bjornstad et al. (2002).
#'
#' @usage data(LBMfromncfpkg)
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
#' clipped so that the region around Japan is central in the dataset. Sampling
#' grid is 0.25 X 0.25 degree.
#'
#' @usage data(Copepod)
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
#' 5, 45–55. \url{https://doi.org/10.5194/essd-5-45-2013}
#'
#' @keywords datasets

"Copepod"
