#' Assessment results for the ICES cod stocks
#'
#' A dataset containing the assessment output for all cod stocks in the ICES areas (except cod-2224) 
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
#' ## load the dataset
#' data(ices.cod)
#' 
#' ## Select the North Sea cod (cod-347d) assessment data.frame
#' nscod <- ices.cod[["cod-347d"]]
#' ## Get the Fmsy from the meta-data
#' fmsy <- as.numeric(attr(nscod, "FMSY"))
#' ## Plot the F/Fmsy of North Sea cod
#' plot(nscod$Year, nscod$F / fmsy, type = "l", xlab = "Year", ylab = "F/Fmsy")
#' lines(nscod$Year, nscod$High_F / fmsy, lty = 2)
#' lines(nscod$Year, nscod$Low_F / fmsy, lty = 2)
#' 
#' ## To see all available meta-data
#' attributes(nscod)
#' @source ICES Stock Assessment Database, 2015/11. ICES, Copenhagen \url{http://http://standardgraphs.ices.dk}
#' @seealso \code{\link{addIces}} for adding ICES assessment results to plots
"ices.cod"

#' Adds ICES assessment results to existing plots
#'
#' Plots lines of asessment results, e.g. F, SSB, etc. to existing plots
#'
#' @param stock The stock code
#' @param col color of the line
#' @param lwd line thickness
#' @param lty line type
#' @param what chr, which column of the result table to plot, see Details.
#' @param mult numeric, multiplier of the data (e.g. 1000 for converting tonnes to kg)
#'
#' @return This function is used for its sidefect. It returns \code{NULL} invisibly
#' @export
#'
#' @note If no active plot exists produces an error
addIces <- function(stock, col="darkgrey", lwd=2, lty = c(2,1,2), what = "FFmsy", mult = 1) {
  ices <- ices.cod[[stock]]
  if(is.null(ices)) stop("Stock ", stock, " was not in the ices.cod dataset")
  nms <- tolower(names(ices))
  what <- tolower(what)
  if(what == "ffmsy") {
    matplot(ices$Year, ices[ , na.omit(pmatch(c("high_f", "f","low_f"), nms, NULL))] / fmsy(ices) * mult, 
            add=TRUE, col=col, lwd = lwd, lty = lty, type="l")
  } else {
    n <- pmatch(what, nms)
    h <- pmatch(paste0("high_", what), nms)
    l <- pmatch(paste0("low_", what), nms)
    matplot(ices$Year, ices[ , c(h, n, l)] * mult, add=TRUE, col=col, lwd = lwd, lty = lty, type="l")
  }
  return(invisible(NULL))
}

fmsy <- function(x) {
  fmsy <- attr(x, "FMSY")
  if(is.null(fmsy)) warning("There is no Fmsy in object x")
  as.numeric(fmsy)
}


