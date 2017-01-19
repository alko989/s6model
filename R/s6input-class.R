#' The s6input class and constructor
#'
#' s6input is an S4 class that contains the input data that are used to 
#' make the assessment. The \code{s6input} function is the constructor of the
#' class. Convienient functions are available to present the data and 
#' \code{plot} them.
#'
#' @slot wf list. 
#' @slot surWF list. 
#' @slot isSimulated logical. 
#' @slot trueParams s6params. 
#' @slot isSurvey logical. 
#' @slot catch numeric. 
#' @slot years numeric. 
#' @slot stockname chacharacter. 
#' @slot stockcode character.
#' @slot speciesname character. 
#'
#' @author alko
#' 
#' @exportClass s6input
#' @name s6input
#' @aliases s6input-class
#' @rdname s6input
#' @include s6params-class.R
#' @export
setClass("s6input",
         slots = c(
           wf = "list",
           surWF = "list", ## this not used at the moment
           isSimulated = "logical",
           trueParams = "s6params",
           isSurvey = "logical",
           catch = "numeric",
           years = "numeric",
           stockname = "character",
           stockcode = "character",
           speciesname = "character"
         ),
         prototype = c(
           wf = list(),
           surWF = list(),
           isSimulated = FALSE,
           trueParams = s6params(),
           isSurvey = FALSE,
           catch = c(),
           years = c(),
           stockname = "generic stock",
           stockcode = "-",
           speciesname = "-"
         )
)


#' @param wf list of data frames with commercial weight frequencies
#' @param surWF list of data frames with scientific survey weight frequencies
#' @param isSimulated logical, are the data simulated?
#' @param trueParams s6params object, the parameters used in the simulation (is used only if isSimulated)
#' @param isSurvey logical, if TRUE the surWF is used and wf is ignored and vice versa
#' @param catch vectors of total catch in tonnes
#' @param years vector of years
#' @param stockname character, the stock name, e.g. North Sea cod
#' @param stockcode character, the stock code, e.g. cod-347d
#' @param speciesname character, the species latin name, e.g. "Gadus morhua"
#'
#' @return A s6input object
#' @export
#' @rdname s6input
#' @include s6params-class.R
s6input <- function(wf = list(),
                    surWF = list(),
                    isSimulated = FALSE,
                    trueParams = s6params(),
                    isSurvey = FALSE,
                    catch = numeric(0),
                    years = numeric(0),
                    stockname = "generic stock",
                    stockcode = "-",
                    speciesname = "-") {
  if (isSurvey) {
    if (length(surWF) != length(catch)) {
      stop("The catch vector has different length from the weight frequency time series")
    }
    if (length(surWF) != length(years)) {
      stop("The years vector has different length from the weight frequency time series")
    }
  } else {
    if (length(wf) != length(catch)) {
      stop("The catch vector has different length from the weight frequency time series")
    }
    if (length(wf) != length(years)) {
      stop("The years vector has different length from the weight frequency time series")
    }
  }
  
  new("s6input", wf = wf, surWF = surWF, isSimulated = isSimulated, 
      trueParams = trueParams, isSurvey = isSurvey, catch = catch, years = years,
      stockname = stockname, stockcode = stockcode, speciesname = speciesname)
}

#' Plot the weight distribution
#'
#' @param x a s6input object 
#' @param ... additional paramters to \code{\link{hist}}
#'
#' @return Nothing, it is just used for to plot the
#' @export
#' @rdname s6input
#' @examples
#' ## simulate a data set
#' d <- simulate(s6params(), binsize = 100, nsim = 100)
#' 
#' ## Plot the data
#' hist(d)
hist.s6input <- function(x, ...) {
  hist(with(x@wf[[1]], rep(Weight, Freq)), main = "", xlab = "Weight", ylab = "Frequency", ...)
}

#' @rdname s6input
#' @export
setMethod('[', signature(x="s6input"), definition=function(x,i){
  res <- x
  if(x@isSurvey) {
    x@surWF <- x@surWF[i]  
  } else {
    x@wf <- x@wf[i]
  }
  x@catch <- x@catch[i]
  x@years <- x@years[i]
  x
})

#' @rdname s6input
#' @export
seq.s6input <- function(x) {
  seq_along(x@years)
}

#' @rdname s6input
#' @export
length.s6input <- function(x) {
  length(x@years)
}

