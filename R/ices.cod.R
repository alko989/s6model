#' Assessment results for the ICES cod stocks
#'
#' A dataset containing the assessment output for all cod stocks in the ICES areas (except cod-2224)#' 
#'
#' @format A list of data frames with columns
#' \itemize{
#'   \item{Year}
#'   \item{Recruitment}
#'   \item{TBiomass}
#'   \item{SSB}
#'   \item{Landings}
#'   \item{Discards}
#'   \item{F}
#' }
#' Meta data are attached to each data frame as attributes. These include the AssessmentYear, the stock name, the stock key and reference points, e.g. Flim, Fpa, FMSY, Bpa
#' @usage data(ices.cod)
#' @details For some of the cod stocks uncertainty bounds of quantity x are given as High_x and Low_x
#' @examples 
#' ## Select the North Sea cod (cod-347d) assessment data.frame
#' nscod <- ices.cod[["cod-347d"]]
#' ## Get the Fmsy from the meta-data
#' fmsy <- as.numeric(attr(nscod, "FMSY"))
#' ## Plot the F/Fmsy of North Sea cod
#' plot(nscod$Year, nscod$F / fmsy, type = "l",)
#' lines(nscod$Year, nscod$High_F / fmsy, lty = 2)
#' lines(nscod$Year, nscod$Low_F / fmsy, lty = 2)
#' 
#' ## To see all available meta-data
#' attributes(nscod)
#' @source ICES Stock Assessment Database, 2015/11. ICES, Copenhagen \url{http://http://standardgraphs.ices.dk}
"ices.cod"