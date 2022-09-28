#' The s6input class constructor
#'
#' s6input is an S3 class that contains the input data that are used to 
#' make the assessment. The \code{s6input} function is the constructor of the
#' class. Convenient functions are available to present the data and 
#' \code{plot} them.
#'
#' @param wf list of data frames with commercial weight frequencies
#' @param surwf list of data frames with scientific survey weight frequencies
#' @param isSimulated logical, are the data simulated?
#' @param trueParams s6params object, the parameters used in the simulation (is used only if isSimulated)
#' @param isSurvey logical, if TRUE the surWF is used and wf is ignored and vice versa
#' @param catch vector of total catch
#' @param catchUnits character, units of catch
#' @param years vector of years
#' @param stockname character, the stock name, e.g. North Sea cod
#' @param stockcode character, the stock code, e.g. cod-347d
#' @param speciesname character, the species latin name, e.g. "Gadus morhua"
#'
#' @author Alexandros Kokkalis <alko@aqua.dtu.dk>
#' 
#' @return A s6input object
#' 
#' @name s6input
#' @aliases s6input-class
#' @rdname s6input
#' @export
#' @examples  
#' ## Empty input
#' inp <- new_s6input()
#' 
s6input <- function(wf = NULL, surwf = NULL, isSimulated = FALSE, trueParameters = NULL,
                    isSurvey = FALSE, catch = NULL, catchUnits = "tonnes", years = NULL, 
                    stockname = "-", stockcode = "-", species = "-") {
  ## Input argument validation
  if (isSurvey) {
    if (length(surwf) != length(catch)) {
      stop("The catch vector has different length from the weight frequency time series")
    }
    if (length(surwf) != length(years)) {
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
  
  e <- new.env()
  e$wf <- wf
  e$surwf <- surwf
  e$isSimulated <- isSimulated
  e$trueParams <- trueParameters
  e$isSurvey <- isSurvey
  e$catch <- catch
  e$catchUnits <- catchUnits
  e$years <- years
  e$stockname <- stockname
  e$stockcode <- stockcode
  e$species <- species
  e$.__env__ <- e
  e$s6version <- s6model:::getVersion("s6model")
  structure(e, class = "s6input")
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
  obs <- with(x$wf[[1]], rep(Weight, Freq))
  hist(obs, main = "", xlab = "Weight", ylab = "Frequency", ...)
}

#' @rdname s6input
#' @export
seq.s6input <- function(x) {
  seq_along(x$years)
}

#' @param x s6input object
#' @rdname s6input
#' @export
length.s6input <- function(x) {
  length(x$years)
}

#' @rdname s6input
#' @export
`[.s6input` <- function(x, i) {
  if (length(x) < max(i)) stop("Out of bounds.", call. = FALSE)
  res <- as.environment(as.list.environment(x, all.names = TRUE))
  res$surwf <- res$surwf[i]  
  res$wf <- res$wf[i]
  res$catch <- res$catch[i]
  res$years <- res$years[i]
  structure(res, class = "s6input")
}

#' @rdname s6input
#' @export
format.s6input <- function(x, ...) {
  args <- list(...)
  objname <- if ("objname" %in% names(args)) args$objname else "inp"
  addline <- function(x, ...) paste0(x, ..., "\n") 
  line <- paste0(rep("-", 80), collapse= "")
  res <- addline("", "Object of class 's6input'")
  res <- addline(res, "Input data for: ", x$stockname, 
                 " (", x$species, ", ", x$stockcode,")")
  res <- addline(res, line)
  years <- x$years
  ny <- length(years)
  yearlines <- ceiling(ny / 10)
  res <- addline(res, "Available years: ", 
                 paste(head(years, 10), collapse = ", "))
  years <- tail(years, -10)
  for(i in seq(yearlines)[-1]) {
    res <- addline(res, "                 ", 
                   #paste(x$years[((i - 1) * 10 + 1):min((i) * 10 , ny)], collapse = ", "))
                   paste(head(years,10), collapse = ", "))
    years <- tail(years, -10)
  }
  res <- addline(res, "Is survey: ", x$isSurvey)
  res <- addline(res, "Is simulated data: ", x$isSurvey)
  res <- addline(res, line)
  res <- addline(res, "Use `plot(", objname, ")` for visualizing the data ",
                 "and `fit(", objname, ")` to make the assessment.")
  res
}

get_res_obj_name <- function(x) {
  find_last <- function(x) {
    if (!is.list(x[[1]]))
      return(x[[1]]) 
    else 
      return(find_last(x[[1]]))
  }
  res <- find_last(x[length(x)])
  if (any(grepl('UseMethod[(]"print"[)]', deparse(x))))
    return("inp")
  if (is.call(res))
    return(deparse(res[[length(res)]]))
}

#' @rdname s6input
#' @export
print.s6input <- function(x, ...) {
  syscall <- sys.calls()
  objname <- get_res_obj_name(syscall)
  cat(format(x, objname = objname), "\n")
}

#' @rdname s6input
#' @export
is.s6input <- function(x) {
  inherits(x, "s6input")
}
