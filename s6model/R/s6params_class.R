#' The s6params class and constructor
#'
#' s6params is an S3 class that contains all model parameters. The 
#' \code{s6params} function is the constructor of the class. Convenient 
#' functions are available to \code{plot}, draw \code{lines}, \code{simulate} 
#' data, return a named list of all parameters (\code{as.list}) and get the 
#' mean parameters of a list of \code{s6params} objects.
#'
#' @param Winf numeric, asymptotic weight
#' @param Fm numeric, fishing mortality
#' @param A numeric, growth parameter
#' @param n numeric, exponent of consumption
#' @param eta_F numeric, 50\% retention size, relative to asymptotic weight
#' @param eta_m numeric, 50\% maturation size, relative to asymptotic weight
#' @param eta_S numeric, 50\% retention size (survey), relative to asymptotic weight
#' @param a numeric, physiological mortality
#' @param epsilon_a numeric, allocation to activity
#' @param epsilon_r numeric, recruitment efficiency
#' @param Wfs numeric, 50\% retention size
#' @param u numeric, selectivity parameter
#' @param ngrid numeric, number of internal weight classes for Fmsy estimation
#' @param wl.a, numeric, weight length relationship multiplier
#' @param wl.b numeric, weight length relationship exponent
#' @param base s6params object
#'
#' @author Alexandros Kokkalis <alko@aqua.dtu.dk>
#' 
#' @name s6params
#' @aliases s6params-class
#' @rdname s6params
#' @export
#' @examples
#' ## Without any arguments gives a \code{s6params} object with default values
#' s6params()
#' 
#' ## An object with some user selected values
#' par1 <- s6params(Winf = 1000, Fm = 0.4, Wfs = 100)
#'  
#' ## Take a s6params object and change one parameter
#' par <- s6params(Winf = (1000), a = (0.4), Fm = (0.2), Wfs = (100))
#' changeMatsize <- par
#' changeMatsize$etam <- (0.3)
#'
#' difference(par, changeMatsize)
#' ##       base comp difference percent.difference
#' ## eta_m 0.25  0.3      -0.05                 20
s6params <- function(
    Winf =      10000,  # Asymptotic weight
    Fm =        0.25,   # Fishing mortality
    A =         4.47,   # Growth parameter
    n =         0.75,   # Exponent of consumption
    eta_m=      0.25,   # Maturation weight (relative to Winf)
    eta_S=      0.001,  # 50% retainment weight, survey gear
    a =         0.22,   # Physiological mortality
    epsilon_a = 0.8,    # Allocation to maintenance
    epsilon_r = 0.1,    # Reproduction efficiency
    Wfs =  500,         # 50% retainment weight of commercial gear
    u =    10,
    ngrid = 1000,       # Grid points for Fmsy estimation
    wl.a = 0.01,        # Weight-length relationship a,b: W = a L^b
    wl.b = 3,                
    base = NULL) {
  ## Get user defined parameters
  mc <- match.call()
  mcl <- as.list(mc[-1])
  userargs <- names(mcl)
  usebase <- "base" %in% userargs
  userargs <- userargs[userargs != "base"]
  
  if(usebase && ! is.s6params(base)) stop("base should be an s6params object", call. = FALSE)
  
  for (n in userargs) {
    if (! is.finite(eval(mcl[[n]]))) stop(paste(n, " should be finite number"), call. = FALSE)
    if (length(eval(mcl[[n]])) != 1) stop(paste(n, " should not be a vector"), call. = FALSE)
  }
  
  e <- new.env()
  e$Winf <- Winf
  e$Fm <- Fm
  e$A <- A
  e$n <- n
  e$eta_m <- eta_m
  e$eta_S <- eta_S
  e$a <- a
  e$epsilon_a <- epsilon_a
  e$epsilon_r <- epsilon_r
  e$Wfs <- Wfs
  e$u <- u
  e$ngrid <- ngrid
  e$wl.a <- wl.a
  e$wl.b <- wl.b
  
  if ("base" %in% names(mcl)) {
    mcl <- mcl[mcl!="base"]
    for (n in setdiff(names(base), userargs)) {
      e[[n]] <- base[[n]]
    }
  }
  
  if (e$Winf <= e$Wfs) {
    warning("The start of fishing occurs at a weight equal or greater than the asymptotic weight")
  }
  
  # e$.__env__ <- e
  e$s6version <- s6model:::getVersion("s6model")
  structure(e, class = "s6params")
}


#' @rdname s6params
#' @export
print.s6params <- function(x, ...) {
  cat(format(x), "\n")
}

#' @rdname s6params
#' @export
is.s6params <- function(x) {
  inherits(x, "s6params")
}

#' @rdname s6params
#' @export
all.equal.s6params <- function(target, current) {
  names <- ls(envir = target)
  res <- TRUE
  for (n in names) {
    if(target[[n]] != current[[n]]) return(FALSE)
  }
  res
}


#' s6params <- function(pars = list()) {
#'   res <- base
#'   nms <- names(pars)
#'   mats <- wfs <- etaf <- 0
#'   for(i in seq(along = nms)) {
#'     nm <- nms[i]
#'     if(nm %in% c("M", "wl.a", "wl.b")) {
#'       slot(res, nm) <- pars[[i]]
#'     } else if (nm == "matSize") {
#'       mats <- i
#'     } else if (nm %in% c("Wfs", "logWfs")) {
#'       wfs <- i
#'     } else if (nm %in% c("eta_F", "logeta_F")) {
#'       etaf <- i
#'     } else if (grepl("^log", nm)) {
#'       slot(res, nm) <- pars[[i]]
#'     } else {
#'       slot(res, paste0("log", nm)) <- log(pars[[i]])
#'     }
#'   }
#'   if(mats > 0) {
#'     res@logeta_m <- log(pars[[mats]]) - res@logWinf
#'     if(any(nms %in% c("logeta_m", "eta_m"))) 
#'       warning("The maturation size is defined both in absolute (matSize) and in relative (eta_m) terms. Only matSize is used.")
#'   }
#'   if (wfs > 0) {
#'     res@logWfs <- if (grepl("^log", nms[wfs])) {
#'       pars[[wfs]]
#'     } else {
#'       log(pars[[wfs]])
#'     }
#'     res@logeta_F <- res@logWfs - res@logWinf
#'     if (etaf > 0)
#'       warning("Do not use Wfs and eta_F at the same time. Only Wfs was used")
#'   }
#'   if (etaf > 0 & wfs == 0) {
#'     res@logeta_F <- if (grepl("^log", nms[etaf])) {
#'       pars[[etaf]]
#'     } else {
#'       log(pars[[etaf]])
#'     }
#'     res@logWfs <- res@logeta_F + res@logWinf
#'   }
#'   if (etaf == 0 & wfs == 0) {
#'     res@logWfs <- res@logeta_F + res@logWinf
#'   }
#'   if (res@logWinf <= res@logWfs) {
#'     warning("The start of fishing occurs at a weight equal or greater than the asymptotic weight")
#'   }
#'   res    
#' }
#' 
#' ##' Calculate the mean value of 
#' ##' @param l list of \code{s6params} objects
#' ##'
#' ##' @return A \code{s6params} object with mean 
#' ##' values of the input \code{s6params}. If the input is NULL it returns NULL
#' ##' and it returns the i
#' ##' 
#' ##' It returns NULL if \code{x} is NULL and 
#' ##' \code{l} if \code{l} is an object of class \code{s6params}
#' ##' @export
#' ##' @rdname s6params
#' meanParameters <- function(l, ...) {
#'   if(is.null(l)) {
#'     warning("Argument l in `meanParameters` is NULL")
#'     return(NULL)
#'   }
#'   if(is(l, "s6params")) return(l)
#'   p <- as.list(s6params())
#'   p$eta_F <- NULL
#'   s6params(setNames(sapply(names(p), function(i) {
#'     allvals <- sapply(l, function(xx) {
#'       if(is.null(xx)) return(NA)
#'       c(as.list(xx)[[i]])  
#'     })
#'     mean(allvals, na.rm = TRUE)
#'   }), names(p)))
#' }
#' 
#' ##' Takes a \code{s6params} object and changes its asymptotic weight
#' ##'
#' ##' The asymptotic weight is changed, along with the relative and absolute sizes of 50\% retention 
#' ##' @param value Numeric. The new asymptotic weight
#' ##' @return \code{s6params} object with changed asymptotic weight, and absolute and
#' ##' relative 50\% retention sizes
#' ##' @author alko
#' ##' @docType methods
#' ##' @rdname Winf-methods
#' ##' @export
#' setGeneric("Winf<-",function(object,value){standardGeneric("Winf<-")})
#' 
#' ##' @rdname Winf-methods
#' ##' @aliases Winf<--methods, Winf<-,s6params-method
#' ##' @name Winfsetter
#' setReplaceMethod(
#'   f = "Winf",
#'   signature = "s6params",
#'   definition = function(object, value) {
#'     object@logWinf <- log(value)
#'     eF <- exp(object@logeta_F)
#'     winf <- exp(object@logWinf)
#'     object@logWfs <- log(eF * winf)
#'     return (object)
#'   })
#' 
#' ##' @param object \code{s6params} object 
#' ##' @export
#' ##' @docType methods
#' ##' @rdname Winf-methods
#' setGeneric("Winf", function(object) standardGeneric("Winf"))
#' 
#' ##' @rdname Winf-methods
#' ##' @aliases Winf,s6params-method
#' setMethod("Winf", 
#'           signature(object = "s6params"), 
#'           function(object) {
#'             exp(object@logWinf)
#'           }
#' )
#' 
#' formatEntry <- function(..., width = 20) {
#'   res <- c()
#'   for(arg in list(...)) {
#'     if(is(arg, "numeric")) {
#'       res <- c(res, round(arg, ifelse(arg < 10, 4, 1)))
#'     } else {
#'       res <- c(res, arg)
#'     }
#'   }
#'   format(paste(res, sep = "", collapse = ""), width = width)
#' }
#' 
#' setMethod("show", "s6params",
#'           function(object) {
#'             width <- min(floor(getOption("width") / 5), 20)
#'             cat(" ___________________________________\n")
#'             cat("|  An object of class 's6params'    |\n")
#'             cat("|___________________________________|",rep("_", width * 3 - 34), "\n", sep="")
#'             cat("|", formatEntry("  Winf  = ", exp(object@logWinf), width = width), 
#'                 "|", formatEntry("  A = ", exp(object@logA), width = width),
#'                 "|", formatEntry("  eps_r = ", exp(object@logepsilon_r), width = width), "|\n",
#'                 "|", formatEntry("  Fm    = ", exp(object@logFm), width = width),
#'                 "|", formatEntry("  a = ", exp(object@loga), width = width),
#'                 "|", formatEntry("  eta_m = ", exp(object@logeta_m), width = width), "|\n",
#'                 "|", formatEntry("  eta_F = ", exp(object@logeta_F), width = width), 
#'                 "|", formatEntry("  n = ", exp(object@logn), width = width),
#'                 "|", formatEntry("  eps_a = ", exp(object@logepsilon_a), width = width),"|\n",
#'                 "|", formatEntry("  eta_S = ", exp(object@logeta_S), width = width),
#'                 "|", formatEntry("  Wfs = ", exp(object@logWfs), width = width),
#'                 "|", formatEntry("  u = ", exp(object@logu), width = width),"|\n",
#'                 "|", formatEntry("  wl.a = ", object@wl.a, width = width),
#'                 "|", formatEntry("  wl.b = ", object@wl.b, width = width), 
#'                 "|", formatEntry("", width = width),"|\n",
#'                 sep="")
#'             cat("|", rep("_", width), "|", rep("_", width), "|", rep("_", width), "|\n", sep="")
#'             cat("\n")
#'           })
#' 
#' 
#' #' @param x a s6params object
#' #' @param xlim the x limits (x1, x2) of the plot.
#' #' @param ... Arguments passed to other methods.
#' #' @note Additional arguments are passed to \code{\link{plot.default}} 
#' #' (from plot) and to \code{\link{lines}} (from lines). From all other 
#' #' functions the extra arguments are ignored.
#' #'
#' #' @export
#' #' @rdname s6params
#' plot.s6params <- function(x, xlim = c(0.001, 1), ...) {
#'   p <- getParams(x)
#'   plot.default(p$w / p$Winf, p$N * (p$w ^ 2), log="xy",
#'                main = "Biomass with respect to relative weight",
#'                xlab = "w/Winf", ylab = "Biomass",
#'                xlim = xlim,
#'                type = "l", ...)
#' }
#' 
#' #' @export
#' #' @rdname s6params
#' lines.s6params <- function(x, ...){
#'   p <- getParams(x)
#'   lines(xy.coords(p$w / p$Winf, p$N * (p$w ^ 2)), ...)
#' }
#' 
#' ##' @export
#' ##' @rdname s6params
#' as.list.s6params <- function(x, ...) {
#'   res <- lapply(slotNames("s6params"), function(nm) {
#'     exp(slot(x, nm))
#'   })
#'   res <- setNames(res, sub(slotNames("s6params"), pattern = "log", replacement = ""))
#'   res$M <- slot(x, "M")
#'   res$wl.a <- slot(x, "wl.a")
#'   res$wl.b <- slot(x, "wl.b")
#'   res
#' }
#' 
#' ##' Difference between two \code{s6params} objects
#' ##' 
#' ##' @param base \code{s6params} object. First object
#' ##' @param comp \code{s6params} object. Second object
#' ##' @return TRUE if they are the same. If there are differences, a data.frame is returned
#' ##' with the untransformed parameter values of the two objects, the relative difference (base - comp)
#' ##' and the percent difference. 
#' ##' @author alko
#' ##' @docType methods
#' ##' @rdname difference-methods
#' ##' @export
#' setGeneric("difference", function(base, comp) {
#'   standardGeneric ("difference")
#' })
#' 
#' ##' @rdname difference-methods
#' ##' @aliases difference,s6params,s6params-method
#' setMethod("difference", c("s6params", "s6params"), function(base, comp) {
#'   res <- sapply(slotNames("s6params"), function(n) {
#'     if (slot(base, n) != slot(comp, n)) {
#'       val1 <- exp(slot(base, n))
#'       val2 <- exp(slot(comp, n))
#'       c(val1, val2, val1 - val2, abs((val1 - val2) / (mean(val1, val2))) * 100)
#'     } else {
#'       c(NA, NA, NA, NA)
#'     }
#'   })
#'   res <- res[,which(!is.na(res[1,])), drop = FALSE]
#'   rownames(res) <- c("base", "comp", "difference", "percent.difference")
#'   colnames(res) <- gsub("log", "", colnames(res))
#'   if(dim(res)[2] == 0) return(FALSE)
#'   round(res, 4)
#' })
#' 
#' ##' Visualizing fit of s6model
#' ##'
#' ##'
#' ##' @title plotFit
#' ##' @name plotFit-methods
#' ##' @aliases plotFit
#' ##' @param object A \code{s6params} object
#' ##' @param data Numeric vector or data.frame with columns Weight and Freq.
#' ##' @param add Boolean. If TRUE, the plot is added to an existing graphics device.
#' ##' @param ... Extra named arguments are passed to the plotting function
#' ##' @return invisible \code{NULL}
#' ##' @docType methods
#' ##' @rdname plotFit-methods
#' ##' @docType methods
#' ##' @export 
#' setGeneric("plotFit", function(object, data, add, ...){ standardGeneric ("plotFit") })
#' 
#' ##' @rdname plotFit-methods
#' ##' @aliases plotFit,s6params,numeric,missing-method
#' setMethod("plotFit", c("s6params", "numeric", "missing"),
#'           function(object, data,...) {plotFit(object, data, FALSE,...)})
#' 
#' ##' @rdname plotFit-methods
#' ##' @aliases plotFit,s6params,data.frame,missing-method
#' setMethod("plotFit", c("s6params", "data.frame", "missing"),
#'           function(object, data,...) {plotFit(object, data, FALSE,...)})
#' 
#' ##' @rdname plotFit-methods
#' ##' @aliases plotFit,s6params,numeric,logical-method
#' setMethod("plotFit", c("s6params", "numeric", "logical"),
#'           function(object, data, add,...) {
#'             p <- getParams(object)
#'             if(add == FALSE) {
#'               plot(p$w, p$pdfN.approx(p$w), type="l", col="blue",
#'                    main="Fitted pdf and histogram of the simulated data",
#'                    xlab="Weight (g)",
#'                    ylab="Probability")
#'             } else {
#'               lines(p$w, p$pdfN.approx(p$w), col="blue", ...)
#'             }
#'             hist(data, freq=FALSE, add=TRUE, breaks="FD")
#'             invisible(NULL)
#'           })
#' 
#' ##' @rdname plotFit-methods
#' ##' @aliases plotFit,s6params,data.frame,logical-method
#' setMethod("plotFit", c("s6params", "data.frame", "logical"),
#'           function(object, data, add, ...) {
#'             p <- getParams(object)
#'             if(add == FALSE) {
#'               plot(p$w, p$pdfN.approx(p$w), type="l",
#'                    main="Fitted pdf and histogram of the simulated data",
#'                    xlab="Weight (g)",
#'                    ylab="Probability")
#'             } else {
#'               lines(p$w, p$pdfN.approx(p$w), col="blue", ...)
#'             }
#'             points(data$Weight, data$Freq/sum(data$Freq)/diff(c(data$Weight,tail(data$Weight,1))),
#'                    pch=".", cex=3)
#'             lines(density(rep(data$Weight, data$Freq)), col=2, lty=2, lwd=2)
#'             ##hist(rep(data$Weight, data$Freq), breaks = 35, add=T, freq=FALSE)
#'             lines(p$w, p$pdfN.approx(p$w), col="blue", lwd = 2, ...)
#'             legend("topright", NULL, c("fitted PDF", "Data kernel density"), col=c("blue","red"),
#'                    lty=1, lwd=2, seg.len=5)
#'             invisible(NULL)
#'           })
#' 
#' ##' Plots growth function
#' ##'
#' ##' @param object A \code{s6params} object
#' ##' @param ... Additional arguments for plot
#' ##' @return Invisible \code{NULL}
#' ##' @author alko
#' ##' @docType methods
#' ##' @rdname plotGrowthMortality
#' ##' @export
#' setGeneric("plotGrowth", function(object, ...) {standardGeneric("plotGrowth")})
#' ##' @aliases plotGrowth,s6params-methods
#' setMethod("plotGrowth", c("s6params"),
#'           function(object, ...) {
#'             p <- getParams(object)
#'             ylim.min <- max(p$g) / 100
#'             ylim.max <- max(p$g)
#'             ylim <- c( ylim.min - ylim.min * 0.01, ylim.max + ylim.max * 0.6)
#'             plot(p$w / p$Winf, p$g, type="n", 
#'                  xlab="",  log="xy", ylab="",
#'                  xlim=c(0.01, 1), ylim = ylim, yaxt="n",xaxt="n", ...)
#'             polygon(c(p$w/p$Winf, 1 ), c(p$psi_m, 0) * ylim.min * 3 + ylim.min, 
#'                     border="lightgrey",col="lightgrey")
#'             lines(p$w / p$Winf, p$g, lwd=3)
#'             abline(v = p$eta_m, lty = 2, lwd = 1.5)
#'             axis(4, at = c(ylim.min, ylim.min * 4), labels = NA, col.axis = "grey28")
#'             mtext(c(0, 100),at=c(ylim.min, ylim.min * 4), side = 4, line = 0.5)
#'             mtext(side=4, at=ylim.min * 2, text="% mature individuals", line = 1.2, col="grey28")
#'             pow <- 1:3
#'             ticksat <- as.vector(sapply(pow, function(p) (2:10) * 10 ^ p))
#'             axis(2, 10 ^ pow, tcl = 0.5, labels = NA)
#'             mtext(10 ^ pow, 2, 0.5, at = 10 ^ pow)
#'             axis(2, ticksat, labels = NA, tcl = 0.25, lwd = 0, lwd.ticks = 1)
#'             title(ylab="Growth rate (g/y)", line = 1.4)
#'             pow <- -2:0
#'             ticksat <- as.vector(sapply(pow, function(p) (2:10) * 10 ^ p))
#'             axis(1, 10^pow, tcl = 0.5, labels = NA)
#'             axis(1, ticksat, labels = NA, tcl = 0.25, lwd = 0, lwd.ticks = 1)
#'             invisible(NULL)
#'           })
#' 
#' ##' Makes a plot of natural and fishing mortalities
#' ##' @return Invisible \code{NULL}
#' ##' @author alko
#' ##' @rdname plotGrowthMortality
#' ##' @docType methods
#' ##' @export
#' setGeneric("plotMortality", function(object, ...){
#'   standardGeneric("plotMortality")
#' })
#' ##' @aliases plotMortality,s6params-method
#' setMethod("plotMortality", c("s6params"),
#'           function(object, ...) {
#'             p <- getParams(object)
#'             plot(p$w / p$Winf, p$m , type = "n", lwd=2, 
#'                  xlab = "",  log = "x", ylab = "", xaxt = "n", yaxt = "n", ...)
#'             title(xlab = expression(w / W[infinity]))
#'             title(ylab = expression(Mortality ~ (y ^ {-1})), line = 1.4)
#'             lines(p$w / p$Winf, p$psi_F * p$Fm, lwd = 3, lty = "dotted")
#'             lines(p$w / p$Winf, p$m - p$psi_F * p$Fm, lty = 2, lwd = 3)
#'             lines(p$w / p$Winf, p$m, lty = 1, lwd = 3)
#'             pow <- -2:0
#'             ticksat <- as.vector(sapply(pow, function(p) (2:10) * 10 ^ p))
#'             axis(1, 10 ^ pow, tcl = 0.5)
#'             axis(1, ticksat, labels = NA, tcl = 0.25, lwd = 0, lwd.ticks = 1)
#'             coords <- par()$usr
#'             ys <- round(seq(0, coords[4], length.out = 4), 2)
#'             axis(2, labels = NA, at = ys, tcl = 0.5)
#'             mtext(ys, side=2, line = 0.5, at = ys)
#'             invisible(NULL)
#'           })
#' 
#' 
#' ##' @param object a \code{s6params} object.
#' ##' @param nsim number of individuals in the simulated sample.
#' ##' @param seed the seed that is passed to \code{\link{set.seed}}.
#' ##' @param binsize numeric, the width of the weight classes in grams.
#' ##' @param keepZeros logical, if TRUE keep bins with zero frequency.
#' ##' @param ndataset integer, the number of datasets to simulate.
#' ##' 
#' ##' @note In \code{simulate}, all zero bins after the last non-zero bin are droped; even if \code{keepZeros} is TRUE.
#' ##'
#' ##' @export
#' ##' @rdname s6params
#' simulate.s6params <- function(object, nsim = 1000, seed = NULL, binsize = 100, 
#'                               keepZeros = TRUE, ndataset = 1, ...) {
#'   if (is.null(seed)) seed <- as.integer(rnorm(1, mean = 1e5, sd = 1e4))
#'   set.seed(seed) 
#'   p <- getParams(object)
#'   l <- seq(binsize / 2, p$Winf - binsize / 2, binsize)
#'   u <- seq(binsize + binsize / 2, p$Winf + binsize / 2, binsize)
#'   pr <- mapply(function(low, up) {
#'     integrate(p$pdfN.approx, lower = low, upper = up)$value},
#'     l, u)
#'   s <- rmultinom(ndataset, nsim, pr)
#'   res <- lapply(seq(ndataset), function(ds) {
#'     df <- data.frame(Weight = seq(binsize, p$Winf + binsize / 2, binsize), 
#'                      Freq = s[, ds])
#'     df <- structure(df, binsize = binsize)
#'     rl <- rle(df$Freq)
#'     df <- head(df, nrow(df) - if(tail(rl$values, 1) == 0) tail(rl$lengths, 1) else 0)
#'     df
#'   })
#'   s6input(wf = res, surWF = list(), isSimulated = TRUE, trueParams = object, 
#'           catch = rep(1, ndataset), ...)
#' }
#' 
#' ##' @export
#' plot.s6input <- function(x, ...) {
#'   plot(x$Weight, x$Freq, type = "h", lwd = 5, ...)
#' }
#' 
#' ##' @export
#' hist.s6input <- function(x, ..., main = "", xlab = "Weight") {
#'   hist(rep(x$Weight, x$Freq), xlab = xlab,  main = main, ...)
#' }
