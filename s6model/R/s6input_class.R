#' The s6input class and constructor
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
#' @param catch vector of total catch
#' @param catchUnits character, units of catch
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
#' ## North Sea cod input data
#' inpall <- s6input(wf = nscoddat, isSurvey = FALSE,
#' stockcode = "cod.27.47d20",
#' species = "Gadus morhua",
#' stockname = "North Sea cod",
#' years = years,
#' catch = catch)
#'
s6input <- function(wf = NULL, surwf = NULL, isSimulated = FALSE, trueParameters = NULL,
                    catch = NULL, catchUnits = "tonnes",
                    stockname = "-", stockcode = "-", species = "-", verbose = options()$verbose) {
  if (is.data.frame(wf) && setequal(names(wf), c("Weight", "Freq"))) {
    wf <- list(wf)
  }
  ## both commercial and survey data
  both <- !is.null(wf) && !is.null(surwf)
  if (both &&
      is.null(names(wf)) && is.null(names(surwf)) && ## no years provided
      length(wf) != length(surwf)) { ## two datasets have different lengths
    stop("Different number of years provided for commercial and survey data, and the years are not specified.")
  }

  if (both && !is.null(names(wf)) && !is.null(surwf)){
    if (!setequal(names(wf), names(surwf))) stop("Not the same years are provided for wf and surwf.")
  }

  if (!is.null(wf) && is.null(names(wf))) {
    if (verbose) message("Years for `wf` not provided, using 1, 2, 3, ... for years.")
    names(wf) <- seq(wf)
  }
  if (!is.null(surwf) && is.null(names(surwf))) {
    if (verbose) message("Years for `surwf` not provided, using 1, 2, 3, ... for years.")
    names(surwf) <- seq(surwf)
  }

  years <- getYears(union(names(wf), names(surwf)))
  nyears <- length(years)

  if (is.null(catch)) {
    if (verbose) message("Catch was not provided.")
    catch <- data.frame(Year = years, Catch = rep(0.012345, nyears))
  }

  stopifnot(is(catch, "data.frame"))
  stopifnot(setequal(names(catch), c("Catch", "Year")))

  catchMissingYears <- setdiff(years, catch$Year)
  if (length(catchMissingYears) != 0)
    stop("Catch was provided, but not for the years: ", paste(catchMissingYears, collapse = ", "))

  e <- new.env()
  e$wf <- wf
  e$surwf <- surwf
  e$isSimulated <- isSimulated
  e$trueParams <- trueParameters
  e$catch <- catch
  e$catchUnits <- catchUnits
  e$years <- years
  e$stockname <- stockname
  e$stockcode <- stockcode
  e$species <- species
  e$.__env__ <- e
  e$s6version <- s6model:::getVersion("s6model")
  structure(e, class = "s6input")

  ## IDEA: lock all s6input environment bindings, make them available with getters and setters
  ## using makeActiveBinding
  ## Example:
  # inp <- s6input()
  #
  # ## Lock binding to catch
  # lockBinding("catch", env = inp)
  # ## Check lock
  # inp$catch  ## works
  # inp$catch <- 12414 ## Fails
  # inp$catch
  #
  # ## Make active binding to access catch
  # makeActiveBinding("getCatch", function(v) {if(missing(v)) cat("get\n") else {cat("set\n"); unlockBinding("catch", env = inp); assign("catch", v, env = inp); lockBinding("catch", env = inp)}; get("catch", env = inp)}, env = inp)
  #
  # inp$getCatch
  # inp$getCatch <- 2 ## Unlocks, changes relocks
  # inp$getCatch
  #
  # inp$catch <- 44 ## Still locked
# inp$catch


}


#' Plot the weight distribution
#'
#' @param x a s6input object
#' @param ... additional paramters to \code{\link{hist}}
#'
#' @export
#' @rdname s6input
#' @examples
#' ## simulate a data set
#' d <- simulate(s6params(), binsize = 100, nsim = 100)
#'
#' ## Plot the data
#' hist(d)
hist.s6input <- function(x, ...) {
  #obs <- with(x$wf[[1]], rep(Weight, Freq))
  obs <- rep(x$wf[[1]][["Weight"]], x$wf[[1]][["Freq"]])
  hist(obs, main = "", xlab = "Weight", ylab = "Frequency", ...)
}

#' @rdname s6input
#' @export
plot.s6input <- function(x, y, ...) {
  hist(x, ...)
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
  res$catch <- res$catch[i, ]
  res$years <- res$years[i]
  structure(res, class = "s6input")
}

#' @rdname s6input
#' @export
format.s6input <- function(x, ...) {
  args <- list(...)
  objname <- if ("objname" %in% names(args)) args$objname else "inp"
  addline <- function(x, ...) paste0(x, ..., "\n")
  line <- paste0(rep("-", 80), collapse = "")
  res <- addline("", "Object of class 's6input'")
  res <- addline(res, "Input data for: ", x$stockname,
                 " (", x$species, ", ", x$stockcode, ")")
  res <- addline(res, line)
  years <- x$years
  ny <- length(years)
  yearlines <- ceiling(ny / 10)
  res <- addline(res, "Available years: ",
                 paste(head(years, 10), collapse = ", "))
  years <- tail(years, -10)
  for (i in seq(yearlines)[-1]) {
    res <- addline(res, "                 ",
                   #paste(x$years[((i - 1) * 10 + 1):min((i) * 10 , ny)], collapse = ", "))
                   paste(head(years, 10), collapse = ", "))
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
#' @importFrom methods is
#' @export
is.s6input <- function(x) {
  inherits(x, "s6input")
}
