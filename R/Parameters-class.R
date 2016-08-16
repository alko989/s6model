#' The Parameters class and constructor
#'
#' Parameters is an S4 class that contains all model parameters. The 
#' \code{parameters} function is a constructor of the class. Convienient 
#' functions are available to \code{plot}, draw \code{lines}, \code{simulate} 
#' data, return a named list of all parameters (\code{as.list}) and get the 
#' mean parameters of a list of \code{Parameter} objects.
#'
#' @section Slots:
#' \describe{
#'   \item{\code{logWinf}:}{Numeric scalar, asymptotic weight}
#'   \item{\code{logFm}:}{Numeric scalar, fishing mortality}
#'   \item{\code{logA}:}{Numeric scalar, growth parameter}
#'   \item{\code{logn}:}{Numeric scalar, exponent of consumption}
#'   \item{\code{logeta_F}:}{Numeric scalar, 50\% retention size, relative to asymptotic weight}
#'   \item{\code{logeta_m}:}{Numeric scalar, 50\% maturation size, relative to asymptotic weight}
#'   \item{\code{logeta_S}:}{Numeric scalar, 50\% retention size (survey), relative to asymptotic weight}
#'   \item{\code{loga}:}{Numeric scalar, physiological mortality}
#'   \item{\code{logepsilon_a}:}{Numeric scalar, allocation to activity}
#'   \item{\code{logepsilon_r}:}{Numeric scalar, recruitment efficiency}
#'   \item{\code{logWfs}:}{Numeric scalar, 50\% retention size}
#'   \item{\code{logu}:}{Numeric scalar, selectivity parameter, width o}
#'   \item{\code{M}:}{Numeric scalar, number of internal weight classes}
#' }
#' @author alko
#' 
#' @exportClass Parameters
#' @name Parameters
#' @aliases Parameters-class
#' @rdname Parameters
#' @export
setClass("Parameters",
         representation(logWinf="numeric",          # Asymptotic weight
                        logFm="numeric",            # Fishing mortality
                        logA="numeric",             # Growth parameter
                        logn="numeric",             # Exponent of consumtion
                        logeta_F="numeric",         # Starting weight of fishing
                                                    # (related to Winf) **Deprecated**
                        logeta_m ="numeric",        # Maturation weight
                                                    # (related to Winf)
                        logeta_S="numeric",         # 
                        loga ="numeric",            # Natural mortality
                        logepsilon_a ="numeric",    # Allocation to maintenance
                        logepsilon_r ="numeric",    # Efficiancy of reproduction
                        logWfs = "numeric",         # Starting weight of fishing
                        logu = "numeric",
                        M = "numeric"),
         prototype(logWinf = log(10000),
                   logFm = log(0.25),
                   logA = log(4.47),
                   logn = log(0.75),
                   logeta_F = log(0.05),
                   logeta_m = log(0.25),
                   logeta_S = log(0.001),
                   loga = log(0.22),
                   logepsilon_a = log(0.8),
                   logepsilon_r = log(0.1),
                   logWfs = log(500),
                   logu=log(10),
                   M = 1000))

#' @param names string vector, parameter names
#' @param vals numeric vector, Values corresponding to the names; log transformed if \code{transformed} is TRUE
#' @param transformed logical, if TRUE vals are the transformed parameter values.
#' @param base \code{Parameters} object, the parameters of this object are used instead of the default values
#' @return The constructor \code{parameters} returns a \code{Parameters} object.
#' @author alko
#' @keywords constructor
#' @examples
#' 
#' ## Without any arguments gives a Parameters object with default values
#' parameters()
#' 
#' ## Changing some parameters gives the corresponding object
#' par1 <- parameters(c("Winf", "Fm", "Wfs"), c(log(1000), log(0.4), log(100)))
#' par2 <- parameters(c("Winf", "Fm", "Wfs"), c(1000 , 0.4, 100), transformed=FALSE)
#'
#' ## Check if the two objects are equal
#' all.equal(par1, par2)
#'
#' ## Take a Parameters object and change one parameter
#' par <- parameters(c("Winf", "a", "Fm", "Wfs"), c(1000, 0.4, 0.2, 100), transformed = FALSE)
#' changeMatsize <- parameters("eta_m", 0.3, transformed = FALSE, base=par)
#'
#' difference(par, changeMatsize)
#' ##       base comp difference percent.difference
#' ## eta_m 0.25  0.3      -0.05                 20
#' @rdname Parameters
#' @export 
parameters <- function(names= c(), vals = c(), transformed = TRUE, base = new("Parameters")) {
  res <- base
  mats <- wfs <- etaf <- 0
  if (length(names) == 1) {
    if (names=="Winf") {
      if (transformed) {
        res@logWinf <- vals
      } else {
        res@logWinf <- log(vals)
      }
      res@logWfs <- res@logeta_F + res@logWinf
      return(res)
    }
  }
  for(i in seq(along=names))
    if(names[i] %in% c("M")) {
      eval(parse(text=paste("res@", names[i]," <- ", vals[i], sep="" )))
    } else if (names[i] == "matSize") {
      mats <- i
    } else if (names[i] == "Wfs") {
      wfs <- i
    } else if (names[i] == "eta_F") {
      etaf <- i
    } else {
      if(transformed) {
        slot(res, paste0("log", names[i])) <- vals[i]
      } else {
        slot(res, paste0("log", names[i])) <- log(vals[i])
      }
    }
  if(mats > 0) {
    res@logeta_m <- log(vals[mats]) - res@logWinf
  }
  if (wfs > 0) {
    if (transformed) {
      res@logWfs <- vals[wfs]
    } else {
      res@logWfs <- log(vals[wfs])
    }
    res@logeta_F <- res@logWfs - res@logWinf
    if (etaf > 0)
      warning("Do not use Wfs and eta_F at the same time. Only Wfs was used")
  }
  if (etaf > 0 & wfs == 0) {
    if (transformed) {
      res@logeta_F <- vals[etaf]
    } else {
      res@logeta_F <- log(vals[etaf])
    }
    res@logWfs <- res@logeta_F + res@logWinf
  }
  if (etaf == 0 & wfs == 0) {
    res@logWfs <- res@logeta_F + res@logWinf
  }
  if (exp(res@logWinf) <= exp(res@logWfs)) {
    warning("The start of fishing occurs at a weight equal or greater than the asymptotic weight")
  }
  res    
}

##' @param x list of Parameters objects
##'
##' @return \code{meanParameters} returns a \code{Parameters} object with mean 
##' values of the input Parameters. It returns NULL if \code{x} is NULL and 
##' \code{x} if \code{x} is an object of class \code{Parameters}
##' @export
##' @rdname Parameters
meanParameters <- function(l) {
  if(is.null(x)) {
    warning("Argument x in `meanParameters` is NULL")
    return(NULL)
  }
  if(is(x, "Parameters")) return(x)
  p <- as.list(parameters())
  do.call(parameters, 
          list(names = names(p), 
               vals = sapply(seq(p), function(i) {
                 allvals <- sapply(x, function(xx) {
                   if(is.null(xx)) return(NA)
                   c(as.list(xx)[[i]])  
                 })
                 mean(allvals, na.rm = TRUE)
               }), transformed = FALSE))
}

##' Takes a Parameters object and changes its asymptotic weight
##'
##' The asymptotic weight is changed, along with the relative and absolute sizes of 50\% retention 
##' @param value Numeric. The new asymptotic weight
##' @return \code{Parameters} object with changed asymptotic weight, and absolute and
##' relative 50\% retention sizes
##' @author alko
##' @docType methods
##' @rdname Winf-methods
##' @export
setGeneric("Winf<-",function(object,value){standardGeneric("Winf<-")})

##' @rdname Winf-methods
##' @aliases Winf<--methods, Winf<-,Parameters-method
##' @name Winfsetter
setReplaceMethod(
  f = "Winf",
  signature = "Parameters",
  definition = function(object, value) {
    object@logWinf <- log(value)
    eF <- exp(object@logeta_F)
    winf <- exp(object@logWinf)
    object@logWfs <- log(eF * winf)
    return (object)
  })

##' @param object \code{Parameters} object 
##' @export
##' @docType methods
##' @rdname Winf-methods
setGeneric("getWinf", function(object) standardGeneric("getWinf"))

##' @rdname Winf-methods
##' @aliases Winf,Parameters-method
setMethod("getWinf", 
          signature(object = "Parameters"), 
          function(object) {
            exp(object@logWinf)
          }
)

formatEntry <- function(..., width = 20) {
  res <- c()
  for(arg in list(...)) {
    if(is(arg, "numeric")) {
      res <- c(res, round(arg, ifelse(arg < 10, 4, 1)))
    } else {
      res <- c(res, arg)
    }
  }
  format(paste(res, sep = "", collapse = ""), width = width)
}

setMethod("show", "Parameters",
          function(object) {
            width <- min(floor(getOption("width") / 5), 20)
            cat(" ___________________________________\n")
            cat("|  An object of class 'Parameters'  |\n")
            cat("|___________________________________|",rep("_", width * 3 - 34), "\n", sep="")
            cat("|", formatEntry("  Winf  = ", exp(object@logWinf), width = width), 
                "|", formatEntry("  A = ", exp(object@logA), width = width),
                "|", formatEntry("  eps_r = ", exp(object@logepsilon_r), width = width), "|\n",
                "|", formatEntry("  Fm    = ", exp(object@logFm), width = width),
                "|", formatEntry("  a = ", exp(object@loga), width = width),
                "|", formatEntry("  eta_m = ", exp(object@logeta_m), width = width), "|\n",
                "|", formatEntry("  eta_F = ", exp(object@logeta_F), width = width), 
                "|", formatEntry("  n = ", exp(object@logn), width = width),
                "|", formatEntry("  eps_a = ", exp(object@logepsilon_a), width = width),"|\n",
                "|", formatEntry("  eta_S = ", exp(object@logeta_S), width = width),
                "|", formatEntry("  Wfs = ", exp(object@logWfs), width = width),
                "|", formatEntry("  u = ", exp(object@logu), width = width),"|\n",
                sep="")
            cat("|", rep("_", width), "|", rep("_", width), "|", rep("_", width), "|\n", sep="")
            cat("\n")
          })


#' @param x a Parameters object
#' @param xlim the x limits (x1, x2) of the plot.
#' @param ... Arguments passed to other methods.
#' @note Additional arguments are passed to \code{\link{plot.default}} (from plot) and to \code{\link{lines}} (from lines). From all other functions the extra arguments are ignored.
#'
#' @export
#' @rdname Parameters
plot.Parameters <- function(x, xlim = c(0.001, 1), ...) {
  p <- getParams(x)
  plot.default(p$w / p$Winf, p$N * (p$w ^ 2), log="xy",
               main="Biomass with respect to relative weight",
               xlab="w/Winf", ylab="Biomass",
               xlim=xlim,
               type="l", ...)
}

#' @export
#' @rdname Parameters
lines.Parameters <- function(x, ...){
  p <- getParams(x)
  lines(xy.coords(p$w / p$Winf, p$N * (p$w ^ 2)), ...)
}

##' @export
##' @rdname Parameters
as.list.Parameters <- function(x, ...) {
  res <- lapply(slotNames("Parameters"), function(nm) {
    exp(slot(x, nm))
  })
  res <- setNames(res, sub(slotNames("Parameters"), pattern = "log", replacement = ""))
  res$M <- slot(x, "M")
  res
}

##' Difference between two \code{Parameters} objects
##' 
##' @param base \code{Parameters} object. First object
##' @param comp \code{Parameters} object. Second object
##' @return TRUE if they are the same. If there are differences, a data.frame is returned
##' with the untransformed parameter values of the two objects, the relative difference (base - comp)
##' and the percent difference 
##' @author alko
##' @docType methods
##' @rdname difference-methods
##' @export
setGeneric("difference", function(base, comp) {
  standardGeneric ("difference")
})

##' @rdname difference-methods
##' @aliases difference,Parameters,Parameters-method
setMethod("difference", c("Parameters", "Parameters"), function(base, comp) {
  res <- data.frame(base=numeric(), comp=numeric(),difference=numeric(), percent.difference=numeric(), stringsAsFactors = FALSE)
  r <- sapply(slotNames("Parameters"), function(n) {
    if (slot(base, n) != slot(comp, n)) {
      val1 <- exp(slot(base, n))
      val2 <- exp(slot(comp, n))
      res[substr(n, 4, nchar(n)), ] <<- c(val1, val2, val1 - val2, abs((val1 - val2) / (mean(val1, val2))) * 100)
    } else {
      NA
    }
  })
  if(dim(res)[1] == 0) return(TRUE)
  round(res, 4)
})

##' Visualizing fit of s6model
##'
##'
##' @title plotFit
##' @name plotFit-methods
##' @aliases plotFit
##' @param object A \code{Parameters} object
##' @param data Numeric vector or data.frame with columns Weight and Freq.
##' @param add Boolean. If TRUE, the plot is added to an existing graphics device.
##' @param ... Extra named arguments are passed to the plotting function
##' @return invisible \code{NULL}
##' @docType methods
##' @rdname plotFit-methods
##' @docType methods
##' @export 
setGeneric("plotFit", function(object, data, add, ...){ standardGeneric ("plotFit") })

##' @rdname plotFit-methods
##' @aliases plotFit,Parameters,numeric,missing-method
setMethod("plotFit", c("Parameters", "numeric", "missing"),
          function(object, data,...) {plotFit(object, data, FALSE,...)})

##' @rdname plotFit-methods
##' @aliases plotFit,Parameters,data.frame,missing-method
setMethod("plotFit", c("Parameters", "data.frame", "missing"),
          function(object, data,...) {plotFit(object, data, FALSE,...)})

##' @rdname plotFit-methods
##' @aliases plotFit,Parameters,numeric,logical-method
setMethod("plotFit", c("Parameters", "numeric", "logical"),
          function(object, data, add,...) {
            p <- getParams(object)
            if(add == FALSE) {
              plot(p$w, p$pdfN.approx(p$w), type="l", col="blue",
                   main="Fitted pdf and histogram of the simulated data",
                   xlab="Weight (g)",
                   ylab="Probability")
            } else {
              lines(p$w, p$pdfN.approx(p$w), col="blue", ...)
            }
            hist(data, freq=FALSE, add=TRUE, breaks="FD")
            invisible(NULL)
          })

##' @rdname plotFit-methods
##' @aliases plotFit,Parameters,data.frame,logical-method
setMethod("plotFit", c("Parameters", "data.frame", "logical"),
          function(object, data, add, ...) {
            p <- getParams(object)
            if(add == FALSE) {
              plot(p$w, p$pdfN.approx(p$w), type="l",
                   main="Fitted pdf and histogram of the simulated data",
                   xlab="Weight (g)",
                   ylab="Probability")
            } else {
              lines(p$w, p$pdfN.approx(p$w), col="blue", ...)
            }
            points(data$Weight, data$Freq/sum(data$Freq)/diff(c(data$Weight,tail(data$Weight,1))),
                   pch=".", cex=3)
            lines(density(rep(data$Weight, data$Freq)), col=2, lty=2, lwd=2)
            ##hist(rep(data$Weight, data$Freq), breaks = 35, add=T, freq=FALSE)
            lines(p$w, p$pdfN.approx(p$w), col="blue", lwd = 2, ...)
            legend("topright", NULL, c("fitted PDF", "Data kernel density"), col=c("blue","red"),
                   lty=1, lwd=2, seg.len=5)
            invisible(NULL)
          })

##' Plots growth function
##'
##' @param object A \code{Parameters} object
##' @param ... Additional arguments for plot
##' @return Invisible \code{NULL}
##' @author alko
##' @docType methods
##' @rdname plotGrowthMortality
##' @export
setGeneric("plotGrowth", function(object, ...) {standardGeneric("plotGrowth")})
##' @aliases plotGrowth,Parameters-methods
setMethod("plotGrowth", c("Parameters"),
          function(object, ...) {
            p <- getParams(object)
            ylim.min <- max(p$g) / 100
            ylim.max <- max(p$g)
            ylim <- c( ylim.min - ylim.min * 0.01, ylim.max + ylim.max * 0.6)
            plot(p$w / p$Winf, p$g, type="n", 
                 xlab="",  log="xy", ylab="",
                 xlim=c(0.01, 1), ylim = ylim, yaxt="n",xaxt="n", ...)
            polygon(c(p$w/p$Winf, 1 ), c(p$psi_m, 0) * ylim.min * 3 + ylim.min, 
                    border="lightgrey",col="lightgrey")
            lines(p$w / p$Winf, p$g, lwd=3)
            abline(v = p$eta_m, lty = 2, lwd = 1.5)
            axis(4, at = c(ylim.min, ylim.min * 4), labels = NA, col.axis = "grey28")
            mtext(c(0, 100),at=c(ylim.min, ylim.min * 4), side = 4, line = 0.5)
            mtext(side=4, at=ylim.min * 2, text="% mature individuals", line = 1.2, col="grey28")
            pow <- 1:3
            ticksat <- as.vector(sapply(pow, function(p) (2:10) * 10 ^ p))
            axis(2, 10 ^ pow, tcl = 0.5, labels = NA)
            mtext(10 ^ pow, 2, 0.5, at = 10 ^ pow)
            axis(2, ticksat, labels = NA, tcl = 0.25, lwd = 0, lwd.ticks = 1)
            title(ylab="Growth rate (g/y)", line = 1.4)
            pow <- -2:0
            ticksat <- as.vector(sapply(pow, function(p) (2:10) * 10 ^ p))
            axis(1, 10^pow, tcl = 0.5, labels = NA)
            axis(1, ticksat, labels = NA, tcl = 0.25, lwd = 0, lwd.ticks = 1)
            invisible(NULL)
          })

##' Makes a plot of natural and fishing mortalities
##' @return Invisible \code{NULL}
##' @author alko
##' @rdname plotGrowthMortality
##' @docType methods
##' @export
setGeneric("plotMortality", function(object, ...){
  standardGeneric("plotMortality")
  })
##' @aliases plotMortality,Parameters-method
setMethod("plotMortality", c("Parameters"),
          function(object, ...) {
            p <- getParams(object)
            plot(p$w / p$Winf, p$m , type = "n", lwd=2, 
                 xlab = "",  log = "x", ylab = "", xaxt = "n", yaxt = "n", ...)
            title(xlab = expression(w / W[infinity]))
            title(ylab = expression(Mortality ~ (y ^ {-1})), line = 1.4)
            lines(p$w / p$Winf, p$psi_F * p$Fm, lwd = 3, lty = "dotted")
            lines(p$w / p$Winf, p$m - p$psi_F * p$Fm, lty = 2, lwd = 3)
            lines(p$w / p$Winf, p$m, lty = 1, lwd = 3)
            pow <- -2:0
            ticksat <- as.vector(sapply(pow, function(p) (2:10) * 10 ^ p))
            axis(1, 10 ^ pow, tcl = 0.5)
            axis(1, ticksat, labels = NA, tcl = 0.25, lwd = 0, lwd.ticks = 1)
            coords <- par()$usr
            ys <- round(seq(0, coords[4], length.out = 4), 2)
            axis(2, labels = NA, at = ys, tcl = 0.5)
            mtext(ys, side=2, line = 0.5, at = ys)
            invisible(NULL)
          })


##' @param object a \code{Parameters} object.
##' @param nsim number of individuals in the simulated sample.
##' @param seed the seed that is passed to \code{\link{set.seed}}.
##' @param binsize numeric, the width of the weight classes in grams.
##' @param keepZeros logical, if TRUE keep bins with zero frequency.
##' @param ndataset integer, the number of datasets to simulate.
##' 
##' @note In \code{simulate}, all zero bins after the last non-zero bin are droped; even if \code{keepZeros} is TRUE.
##'
##' @export
##' @rdname Parameters
simulate.Parameters <- function(object, nsim = 1000, seed = NULL, binsize = 100, 
                                keepZeros = TRUE, ndataset = 1, ...) {
  if ( ! is.null(seed)) set.seed(seed)
  p <- getParams(object)
  l <- seq(binsize / 2, p$Winf - binsize / 2, binsize)
  u <- seq(binsize + binsize / 2, p$Winf + binsize / 2, binsize)
  pr <- mapply(function(low, up) {
    integrate(p$pdfN.approx, lower = low, upper = up)$value},
    l, u)
  s <- rmultinom(ndataset, nsim, pr)
  res <- lapply(seq(ndataset), function(ds) {
    df <- data.frame(Weight = seq(binsize, p$Winf + binsize / 2, binsize), 
                     Freq = s[, ds])
    df <- structure(df, binsize = binsize)
    rl <- rle(df$Freq)
    df <- head(df, nrow(df) - if(tail(rl$values, 1) == 0) tail(rl$lengths, 1) else 0)
    class(df) <-c("s6Input", "data.frame")  
    df
  })
  if(ndataset == 1) return(res[[1]])
  res
}

##' @export
plot.s6Input <- function(x, ...) {
  plot(x, ...)
}

##' @export
hist.WeightFreq <- function(x, ..., main = "", xlab = "Weight") {
  hist(rep(x$Weight, x$Freq), xlab = xlab,  main = main, ...)
}