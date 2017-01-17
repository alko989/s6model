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
#' @export s6input
s6input <- setClass("s6input",
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
                      comWF = list(),
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



