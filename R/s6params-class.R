#' The s6params class and constructor
#'
#' s6params is an S4 class that contains all model parameters. The 
#' \code{s6params} function is a constructor of the class. Convienient 
#' functions are available to \code{plot}, draw \code{lines}, \code{simulate} 
#' data, return a named list of all parameters (\code{as.list}) and get the 
#' mean parameters of a list of \code{s6params} objects.
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
#' @exportClass s6params
#' @name s6params
#' @aliases s6params-class
#' @rdname s6params
#' @export
setClass("s6params",
         representation(
           logWinf="numeric",          # Asymptotic weight
           logFm="numeric",            # Fishing mortality
           logA="numeric",             # Growth parameter
           logn="numeric",             # Exponent of consumtion
           logeta_F="numeric",         # 50% ratainment weight, commercial gear
           # (related to Winf) **Deprecated**
           logeta_m ="numeric",        # Maturation weight (relative to Winf)
           logeta_S="numeric",         # 50% retainment weight, survey gear
           loga ="numeric",            # Natural mortality
           logepsilon_a ="numeric",    # Allocation to maintenance
           logepsilon_r ="numeric",    # Efficiancy of reproduction
           logWfs = "numeric",         # Starting weight of fishing
           logu = "numeric",
           M = "numeric"),
         prototype(
           logWinf = log(10000),
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
           M = 1000)
         )

#' @param pars named list, parameter values. If the name starts with log the values are expected to be log transformed
#' @param base \code{s6params} object, the parameters of this object are used instead of the default values
#' @return A \code{s6params} object.
#' @author alko
#' @keywords constructor
#' @examples
#' 
#' ## Without any arguments gives a \code{s6params} object with default values
#' s6params()
#' 
#' ## Log transformed s6params are expected if the parameter name starts with 'log'
#' par1 <- s6params(c(logWinf = log(1000), logFm = log(0.4), logWfs = log(100)))
#' 
#' ## The same is achieved with not transformed parameters
#' par2 <- s6params(c(Winf = 1000, Fm = 0.4, Wfs = 100))
#'
#' ## Take a s6params object and change one parameter
#' par <- s6params(list(Winf = 1000, a = 0.4, Fm = 0.2, Wfs = 100))
#' changeMatsize <- s6params(list(eta_m = 0.3), base= par)
#'
#' difference(par, changeMatsize)
#' ##       base comp difference percent.difference
#' ## eta_m 0.25  0.3      -0.05                 20
#' @rdname s6params
#' @export 
s6params <- function(pars = list(), base = new("s6params")) {
  res <- base
  nms <- names(pars)
  mats <- wfs <- etaf <- 0
  for(i in seq(along = nms)) {
    nm <- nms[i]
    if(nm %in% c("M")) {
      slot(res, nm) <- pars[[i]]
    } else if (nm == "matSize") {
      mats <- i
    } else if (nm %in% c("Wfs", "logWfs")) {
      wfs <- i
    } else if (nm %in% c("eta_F", "logeta_F")) {
      etaf <- i
    } else if (grepl("^log", nm)) {
      slot(res, nm) <- pars[[i]]
    } else {
      slot(res, paste0("log", nm)) <- log(pars[[i]])
    }
  }
  if(mats > 0) {
    res@logeta_m <- log(vals[mats]) - res@logWinf
  }
  if (wfs > 0) {
    res@logWfs <- if (grepl("^log", nms[wfs])) {
      pars[[wfs]]
    } else {
      log(pars[[wfs]])
    }
    res@logeta_F <- res@logWfs - res@logWinf
    if (etaf > 0)
      warning("Do not use Wfs and eta_F at the same time. Only Wfs was used")
  }
  if (etaf > 0 & wfs == 0) {
    res@logeta_F <- if (grepl("^log", nms[etaf])) {
      pars[[etaf]]
    } else {
      log(pars[[etaf]])
    }
    res@logWfs <- res@logeta_F + res@logWinf
  }
  if (etaf == 0 & wfs == 0) {
    res@logWfs <- res@logeta_F + res@logWinf
  }
  if (res@logWinf <= res@logWfs) {
    warning("The start of fishing occurs at a weight equal or greater than the asymptotic weight")
  }
  res    
}

##' Calculate the mean value of 
##' @param l list of \code{s6params} objects
##'
##' @return A \code{s6params} object with mean 
##' values of the input \code{s6params}. If the input is NULL it returns NULL
##' and it returns the i
##' 
##' It returns NULL if \code{x} is NULL and 
##' \code{l} if \code{l} is an object of class \code{s6params}
##' @export
##' @rdname s6params
meanParameters <- function(l, ...) {
  if(is.null(l)) {
    warning("Argument l in `meanParameters` is NULL")
    return(NULL)
  }
  if(is(l, "s6params")) return(l)
  p <- as.list(s6params())
  do.call(s6params, 
          list(names = names(p), 
               vals = sapply(seq(p), function(i) {
                 allvals <- sapply(l, function(xx) {
                   if(is.null(xx)) return(NA)
                   c(as.list(xx)[[i]])  
                 })
                 mean(allvals, na.rm = TRUE)
               }), transformed = FALSE))
}

##' Takes a \code{s6params} object and changes its asymptotic weight
##'
##' The asymptotic weight is changed, along with the relative and absolute sizes of 50\% retention 
##' @param value Numeric. The new asymptotic weight
##' @return \code{s6params} object with changed asymptotic weight, and absolute and
##' relative 50\% retention sizes
##' @author alko
##' @docType methods
##' @rdname Winf-methods
##' @export
setGeneric("Winf<-",function(object,value){standardGeneric("Winf<-")})

##' @rdname Winf-methods
##' @aliases Winf<--methods, Winf<-,s6params-method
##' @name Winfsetter
setReplaceMethod(
  f = "Winf",
  signature = "s6params",
  definition = function(object, value) {
    object@logWinf <- log(value)
    eF <- exp(object@logeta_F)
    winf <- exp(object@logWinf)
    object@logWfs <- log(eF * winf)
    return (object)
  })

##' @param object \code{s6params} object 
##' @export
##' @docType methods
##' @rdname Winf-methods
setGeneric("Winf", function(object) standardGeneric("Winf"))

##' @rdname Winf-methods
##' @aliases Winf,s6params-method
setMethod("Winf", 
          signature(object = "s6params"), 
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

setMethod("show", "s6params",
          function(object) {
            width <- min(floor(getOption("width") / 5), 20)
            cat(" ___________________________________\n")
            cat("|  An object of class 's6params'    |\n")
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


#' @param l a s6params object
#' @param xlim the x limits (x1, x2) of the plot.
#' @param ... Arguments passed to other methods.
#' @note Additional arguments are passed to \code{\link{plot.default}} 
#' (from plot) and to \code{\link{lines}} (from lines). From all other 
#' functions the extra arguments are ignored.
#'
#' @export
#' @rdname s6params
plot.s6params <- function(x, xlim = c(0.001, 1), ...) {
  p <- getParams(x)
  plot.default(p$w / p$Winf, p$N * (p$w ^ 2), log="xy",
               main = "Biomass with respect to relative weight",
               xlab = "w/Winf", ylab = "Biomass",
               xlim = xlim,
               type = "l", ...)
}

#' @export
#' @rdname s6params
lines.s6params <- function(x, ...){
  p <- getParams(x)
  lines(xy.coords(p$w / p$Winf, p$N * (p$w ^ 2)), ...)
}

##' @export
##' @rdname s6params
as.list.s6params <- function(x, ...) {
  res <- lapply(slotNames("s6params"), function(nm) {
    exp(slot(x, nm))
  })
  res <- setNames(res, sub(slotNames("s6params"), pattern = "log", replacement = ""))
  res$M <- slot(x, "M")
  res
}

##' Difference between two \code{s6params} objects
##' 
##' @param base \code{s6params} object. First object
##' @param comp \code{s6params} object. Second object
##' @return TRUE if they are the same. If there are differences, a data.frame is returned
##' with the untransformed parameter values of the two objects, the relative difference (base - comp)
##' and the percent difference. 
##' @author alko
##' @docType methods
##' @rdname difference-methods
##' @export
setGeneric("difference", function(base, comp) {
  standardGeneric ("difference")
})

##' @rdname difference-methods
##' @aliases difference,s6params,s6params-method
setMethod("difference", c("s6params", "s6params"), function(base, comp) {
  res <- sapply(slotNames("s6params"), function(n) {
    if (slot(base, n) != slot(comp, n)) {
      val1 <- exp(slot(base, n))
      val2 <- exp(slot(comp, n))
       c(val1, val2, val1 - val2, abs((val1 - val2) / (mean(val1, val2))) * 100)
    } else {
      c(NA, NA, NA, NA)
    }
  })
  res <- res[,which(!is.na(res[1,])), drop = FALSE]
  rownames(res) <- c("base", "comp", "difference", "percent.difference")
  colnames(res) <- gsub("log", "", colnames(res))
  if(dim(res)[2] == 0) return(TRUE)
  round(res, 4)
})

##' Visualizing fit of s6model
##'
##'
##' @title plotFit
##' @name plotFit-methods
##' @aliases plotFit
##' @param object A \code{s6params} object
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
##' @aliases plotFit,s6params,numeric,missing-method
setMethod("plotFit", c("s6params", "numeric", "missing"),
          function(object, data,...) {plotFit(object, data, FALSE,...)})

##' @rdname plotFit-methods
##' @aliases plotFit,s6params,data.frame,missing-method
setMethod("plotFit", c("s6params", "data.frame", "missing"),
          function(object, data,...) {plotFit(object, data, FALSE,...)})

##' @rdname plotFit-methods
##' @aliases plotFit,s6params,numeric,logical-method
setMethod("plotFit", c("s6params", "numeric", "logical"),
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
##' @aliases plotFit,s6params,data.frame,logical-method
setMethod("plotFit", c("s6params", "data.frame", "logical"),
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
##' @param object A \code{s6params} object
##' @param ... Additional arguments for plot
##' @return Invisible \code{NULL}
##' @author alko
##' @docType methods
##' @rdname plotGrowthMortality
##' @export
setGeneric("plotGrowth", function(object, ...) {standardGeneric("plotGrowth")})
##' @aliases plotGrowth,s6params-methods
setMethod("plotGrowth", c("s6params"),
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
##' @aliases plotMortality,s6params-method
setMethod("plotMortality", c("s6params"),
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


##' @param object a \code{s6params} object.
##' @param nsim number of individuals in the simulated sample.
##' @param seed the seed that is passed to \code{\link{set.seed}}.
##' @param binsize numeric, the width of the weight classes in grams.
##' @param keepZeros logical, if TRUE keep bins with zero frequency.
##' @param ndataset integer, the number of datasets to simulate.
##' 
##' @note In \code{simulate}, all zero bins after the last non-zero bin are droped; even if \code{keepZeros} is TRUE.
##'
##' @export
##' @rdname s6params
simulate.s6params <- function(object, nsim = 1000, seed = NULL, binsize = 100, 
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
    class(df) <-c("s6input", "data.frame")  
    df
  })
  if(ndataset == 1) return(res[[1]])
  res
}

##' @export
plot.s6input <- function(x, ...) {
  plot(x$Weight, x$Freq, type = "h", lwd = 5, ...)
}

##' @export
hist.s6input <- function(x, ..., main = "", xlab = "Weight") {
  hist(rep(x$Weight, x$Freq), xlab = xlab,  main = main, ...)
}