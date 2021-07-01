#' The Parameters class and constructor
#'
#' A class that contains all model parameters
#'
#' @section Slots:
#' \describe{
#'   \item{\code{logWinf}:}{Numeric scalar. Asymptotic weight}
#'   \item{\code{logFm}:}{Numeric scalar. Fishing mortality}
#'   \item{\code{logA}:}{Numeric scalar. Growth parameter}
#'   \item{\code{logn}:}{Numeric scalar. Exponent of consumption}
#'   \item{\code{logeta_F}:}{Numeric scalar. 50\% retention size, relative to asymptotic weight}
#'   \item{\code{logeta_m}:}{Numeric scalar. 50\% maturation size, relative to asymptotic weight}
#'   \item{\code{logeta_S}:}{Numeric scalar. 50\% retention size (survey), relative to asymptotic weight}
#'   \item{\code{loga}:}{Numeric scalar. Physiological mortality}
#'   \item{\code{logepsilon_a}:}{Numeric scalar. Allocation to activity}
#'   \item{\code{logepsilon_r}:}{Numeric scalar. Recruitment efficiency}
#'   \item{\code{logWfs}:}{Numeric scalar. 50\% retention size}
#'   \item{\code{logu}:}{Numeric scalar. Selectivity parameter, width o}
#'   \item{\code{M}:}{Numeric scalar. Number of internal weight classes}
#'   \item{\code{scaleWinf}:}{Numeric scalar. Scale of asymptotic weight}
#'   \item{\code{scaleFm}:}{Numeric scalar. Scale of fishing mortality}
#'   \item{\code{scaleA}:}{Numeric scalar. Scale of growth parameter}
#'   \item{\code{scalen}:}{Numeric scalar. Scale of exponent of consumption}
#'   \item{\code{scaleeta_F}:}{Numeric scalar. Scale of 50\% retantion size}
#'   \item{\code{scaleeta_m}:}{Numeric scalar. Scale of 50\% maturation size}
#'   \item{\code{scaleeta_S}:}{Numeric scalar. Scale of survey gear 50\% retantion size}
#'   \item{\code{scalea}:}{Numeric scalar. Scale of the physiological mortality}
#'   \item{\code{scaleepsilon_a}:}{Numeric scalar. Scale of the allocation to activity}
#'   \item{\code{scaleepsilon_r}:}{Numeric scalar. Scale of recruitment efficiency}
#'   \item{\code{scaleWfs}:}{Numeric scalar. Scale of size of 50\% retention}
#'   \item{\code{scaleu}:}{Numeric scalar.}
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
                        logeta_S="numeric",
                        loga ="numeric",            # Natural mortality
                        logepsilon_a ="numeric",    # Allocation to maintenance
                        logepsilon_r ="numeric",    # Efficiancy of reproduction
                        logWfs = "numeric",         # Starting weight of fishing
                        logu = "numeric",           # Width of change from 0 to F fishing mortality
                        M ="numeric",
                        scaleWinf="numeric",        # Asymptotic weight
                        scaleFm="numeric",          # Fishing mortality
                        scaleA="numeric",           # Growth parameter
                        scalen="numeric",           # Exponent of consumtion
                        scaleeta_F="numeric",       # Starting weight of fishing
                                                    #(related to Winf)
                        scaleeta_m ="numeric",      # Maturation weight
                                                    # (related to Winf)
                        scalea ="numeric",          # Natural mortality
                        scaleepsilon_a ="numeric",  # Allocation to maintenance
                        scaleepsilon_r ="numeric",
                        scaleeta_S="numeric",
                        scaleWfs = "numeric",
                        scaleu = "numeric"),
         prototype(logWinf = log(2),
                   logFm = log(1),
                   logA = log(1),
                   logn = log(1),
                   logeta_F = log(1),
                   logeta_m = log(1),
                   logeta_S = log(1),
                   loga = log(1),
                   logepsilon_a = log(1),
                   logepsilon_r = log(1),
                   logWfs = log(1),
                   logu=log(1),
                   M = 1000,
                   scaleWinf = 10000,
                   scaleFm = 0.25,
                   scaleA = 4.47,
                   scalen = 0.75,
                   scaleeta_F = 0.05,
                   scaleeta_m = 0.25,
                   scalea = 0.22,
                   scaleepsilon_a = 0.8,
                   scaleepsilon_r = 0.1,
                   scaleeta_S=0.001,
                   scaleWfs = 1000,
                   scaleu=10))

setGeneric("getscaleWinf",function(object){standardGeneric ("getscaleWinf")})
setMethod("getscaleWinf","Parameters", function(object){ return(object@scaleWinf) })
setGeneric("getscaleFm",function(object){standardGeneric ("getscaleFm")})
setMethod("getscaleFm","Parameters", function(object){ return(object@scaleFm) })
setGeneric("getscaleA",function(object){standardGeneric ("getscaleA")})
setMethod("getscaleA","Parameters", function(object){ return(object@scaleA) })
setGeneric("getscalen",function(object){standardGeneric ("getscalen")})
setMethod("getscalen","Parameters", function(object){ return(object@scalen) })
setGeneric("getscaleeta_F",function(object){standardGeneric ("getscaleeta_F")})
setMethod("getscaleeta_F","Parameters", function(object){ return(object@scaleeta_F) })
setGeneric("getscaleeta_m",function(object){standardGeneric ("getscaleeta_m")})
setMethod("getscaleeta_m","Parameters", function(object){ return(object@scaleeta_m) })
setGeneric("getscalea",function(object){standardGeneric ("getscalea")})
setMethod("getscalea","Parameters", function(object){ return(object@scalea) })
setGeneric("getscaleepsilon_a",function(object){standardGeneric ("getscaleepsilon_a")})
setMethod("getscaleepsilon_a","Parameters", function(object){ return(object@scaleepsilon_a) })
setGeneric("getscaleepsilon_r",function(object){standardGeneric ("getscaleepsilon_r")})
setMethod("getscaleepsilon_r","Parameters", function(object){ return(object@scaleepsilon_r) })
setGeneric("getscaleeta_S",function(object){standardGeneric ("getscaleeta_S")})
setMethod("getscaleeta_S","Parameters", function(object){ return(object@scaleeta_S) })
setGeneric("getscaleWfs",function(object){standardGeneric ("getscaleWfs")})
setMethod("getscaleWfs","Parameters", function(object){ return(object@scaleWfs) })
setGeneric("getscaleu",function(object){standardGeneric ("getscaleu")})
setMethod("getscaleu","Parameters", function(object){ return(object@scaleu) })

#' @param names String vector. Contains the names of the parameters that will
#' have non default values.
#' @param vals Numeric vector. The corresponding values, transformed if \code{transformed} is TRUE.
#' @param transformed Boolean. If TRUE vals should contain the transformed parameter values.
#' @param base \code{Parameters} object. The parameter values will be used instead of the default values.
#' @return Returns an object of the Parameters class
#' @author alko
#' @keywords constructor
#' @examples
#' 
#' ## Without any arguments gives a Parameters object with default values
#' parameters()
#' 
#' ## Changing some parameters gives the corresponding object
#' par1 <- parameters(c("Winf", "Fm", "Wfs"), c(log(1000 / 10000), log(0.4 / 0.25), log(100 / 1000)))
#' par2 <- parameters(c("Winf", "Fm", "Wfs"), c(1000 , 0.4, 100), transformed=FALSE)
#'
#' ## Check if the two objects are equal
#' all.equal(par1, par2)
#'
#' ## Take a Parameters object and change one parameter
#' par <- parameters(c("Winf", "a", "Fm", "Wfs"), c(1000, 0.4, 0.2, 100), transformed = FALSE)
#' changeMatsize <- parameters("eta_m", 0.3, transformed =FALSE, base=par)
#'
#' difference(par, changeMatsize)
#' ##       base comp difference percent.difference
#' ## eta_m 0.25  0.3      -0.05                 20
#' @rdname Parameters
#' @export 
parameters <- function(names= c(), vals = c(), transformed=TRUE, base=new("Parameters"))
{
  res <- base
  mats <- wfs <- etaf <- 0
  if(length(names) == 1)  {
    if(names=="Winf") {
      if(transformed) {
        res@logWinf <- vals
      } else {
        res@logWinf <- log(vals / getscaleWinf(res))
      }
      res@logWfs <- log(exp(res@logeta_F) * res@scaleeta_F * exp(res@logWinf) * res@scaleWinf/res@scaleWfs)
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
        eval(parse(text=paste("res@log", names[i]," <- ", vals[i] , sep="" )))
      } else {
        scale <- do.call(paste0("getscale", names[i]), list(res))
        eval(parse(text=paste0("res@log", names[i], " <- " , log(vals[i] / scale))))
      }
    }
  if(mats > 0) {
    res@logeta_m <- log(vals[mats] / (exp(res@logWinf)*res@scaleWinf) / res@scaleeta_m)
  }
  if(wfs > 0) {
    if(transformed) {
      res@logWfs <- vals[wfs]
    } else {
      res@logWfs <- log(vals[wfs] /getscaleWfs(res))
    }
    res@logeta_F <- log((exp(res@logWfs) * res@scaleWfs) / (exp(res@logWinf) * res@scaleWinf) / res@scaleeta_F)
    if(etaf > 0)
      warning("Do not use Wfs and eta_F at the same time. Only Wfs was used")
  }
  if(etaf>0 & wfs==0) {
    if(transformed) {
      res@logeta_F <-  vals[etaf]
    } else {
      res@logeta_F <- log(vals[etaf] / getscaleeta_F(res))
    }
    res@logWfs <- log(exp(res@logeta_F) * res@scaleeta_F* exp(res@logWinf) * res@scaleWinf/res@scaleWfs)
  }
  if(etaf == 0 & wfs ==0) {
    res@logWfs <- log(exp(res@logeta_F) * res@scaleeta_F * exp(res@logWinf) * res@scaleWinf/res@scaleWfs)
  }
  if(exp(res@logWinf)*res@scaleWinf <= exp(res@logWfs)*res@scaleWfs)
    warning("The start of fishing occurs at a weight equal or greater than the asymptotic weight")
  res    
}

##' @export
##' @rdname Parameters
meanParameters <- function(x) {
  if(is.null(x)) {
    warning("Argument x in `meanParameters` is NULL")
    return(NULL)
  }
  p <- as.list(parameters())
  do.call(parameters, 
          list(names = names(p), 
               vals = sapply(seq(p), function(i) mean(sapply(x, function(xx) {
                 if(is.null(xx)) return(NA)
                 c(as.list(xx)[[i]])  
               }), na.rm = TRUE)), transformed = FALSE))
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
  f="Winf",
  signature="Parameters",
  definition=function(object,value){
    object@logWinf <- log(value/object@scaleWinf)
    eF <- exp(object@logeta_F) * object@scaleeta_F
    winf <- exp(object@logWinf)* object@scaleWinf
    object@logWfs <- log(eF*winf/object@scaleWfs)
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
            exp(object@logWinf) * object@scaleWinf
          }
)

formatEntry <- function(..., width=20) {
  res <- c()
  for(arg in list(...)) {
    if(is(arg, "numeric")) {
      res <- c(res,round(arg, ifelse(arg < 10, 4, 1)))
    } else {
      res <- c(res, arg)
    }
  }
  format(paste(res, sep="", collapse=""), width=width)
}
setMethod("show", "Parameters",
          function(object) {
            width <- min(floor(getOption("width") / 5), 20)
            cat(" ___________________________________\n")
            cat("|  An object of class 'Parameters'  |\n")
            cat("|___________________________________|",rep("_", width * 3 - 34), "\n", sep="")
            cat("|", formatEntry("  Winf  = ",exp(object@logWinf)*object@scaleWinf, width = width), 
                "|", formatEntry("  A = ", exp(object@logA) * object@scaleA, width = width),
                "|", formatEntry("  eps_r = ",exp(object@logepsilon_r)*object@scaleepsilon_r, width = width), "|\n",
                "|", formatEntry("  Fm    = ", exp(object@logFm) * object@scaleFm, width = width),
                "|", formatEntry("  a = ", exp(object@loga) * object@scalea, width = width),
                "|", formatEntry("  eta_m = ", exp(object@logeta_m) * object@scaleeta_m, width = width), "|\n",
                "|", formatEntry("  eta_F = ", exp(object@logeta_F) * object@scaleeta_F, width = width), 
                "|", formatEntry("  n = ", exp(object@logn) * object@scalen, width = width),
                "|", formatEntry("  eps_a = ",   exp(object@logepsilon_a) * object@scaleepsilon_a, width = width),"|\n",
                "|", formatEntry("  eta_S = ",exp(object@logeta_S)*object@scaleeta_S, width = width),
                "|", formatEntry("  Wfs = ", exp(object@logWfs) * object@scaleWfs, width = width),
                "|", formatEntry("  u = ", exp(object@logu) * object@scaleu, width = width),"|\n",
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
plot.Parameters <- function(x, xlim=c(0.001, 1), ...) {
  p <- getParams(x)
  plot.default(p$w / p$Winf, p$N*(p$w^2), log="xy",
               main="Biomass with respect to relative weight",
               xlab="w/Winf", ylab="Biomass",
               xlim=xlim,
               type="l", ...)
}

#' @export
#' @rdname Parameters
lines.Parameters <- function(x, ...){
  p <- getParams(x)
  lines(xy.coords(p$w / p$Winf, p$N*(p$w^2)), ...)
}

##' @export
##' @rdname Parameters
as.list.Parameters <- function(x, ...) {
  with(getParams(x), list(Winf=Winf, Fm=Fm, Wfs=Wfs, eta_m=eta_m, epsilon_r=epsilon_r, epsilon_a=epsilon_a, A=A, a=a, n=n, u=u))
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
  res <- data.frame(base=numeric(), comp=numeric(),difference=numeric(),percent.difference=numeric(), stringsAsFactors=FALSE)
  r <- sapply(slotNames("Parameters"), function(n) {
    if(eval(parse(text=paste("base@" , n, sep=""))) !=
         eval(parse(text=paste("comp@" , n, sep="")))) {
      val1 <- exp(eval(parse(text=paste("base@" , n, sep="")))) * eval(parse(text=paste("base@scale" , substr(n, 4, nchar(n)), sep="")))
      val2 <- exp(eval(parse(text=paste("comp@" , n, sep="")))) * eval(parse(text=paste("comp@scale" , substr(n, 4, nchar(n)), sep="")))
      res[substr(n, 4, nchar(n)), ] <<- c(val1,val2, val1 - val2,abs((val1-val2)/(mean(val1,val2)))*100)
    } else {
      NA
    }
  })
  if(dim(res)[1] == 0) return( TRUE)
  round(res,4)
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
              plot(p$w, p$pdfN.approx(p$w), type="l", 
                   # main="Fitted pdf and histogram of the simulated data",col="blue",
                   xlab="Weight (g)",
                   ylab="Probability")
            } else {
              lines(p$w, p$pdfN.approx(p$w), col="grey28", ...)
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
                   # main="Fitted pdf and histogram of the simulated data",
                   xlab="Weight (g)",
                   ylab="Probability")
            } else {
              lines(p$w, p$pdfN.approx(p$w), col="grey28", ...)
            }
            points(data$Weight, data$Freq/sum(data$Freq)/diff(c(data$Weight,tail(data$Weight,1))),
                   pch=16, cex=1, col="#7E6148B2") #
            # lines(density(rep(data$Weight, data$Freq)), col=2, lty=2, lwd=2)
            ##hist(rep(data$Weight, data$Freq), breaks = 35, add=T, freq=FALSE)
            lines(p$w, p$pdfN.approx(p$w), lwd = 2, ...) #col="blue", 
            legend("topright", NULL, c("Fitted Probability Density Function (PDF)"),
                   lty=1, lwd=2, seg.len=5, bty = "n") #, "Data kernel density" ,"red" col=c("blue"),
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
            # Set y axis from 1 to 100
            ylim.min <- max(p$g) / 100
            ylim.max <- max(p$g)
            #  plot the growth function. x limits are 0.1 to 1, y limits are a little extended
            plot(p$w / p$Winf, p$g, type="n", 
                 xlab="",  log="xy", ylab="",
                 xlim=c(0.01, 1), ylim=c( ylim.min- ylim.min*0.01, ylim.max + ylim.max*0.6), yaxt="n",xaxt="n", ...)
            # Plot the selectivity curve as a polygon. 
            polygon(c(p$w/p$Winf, 1 ), c(p$psi_m, 0) * ylim.min * 100 + ylim.min, border="lightgrey",col="lightgrey") #ylim.min * 3 + ylim.min
            title(xlab=expression(w/W[infinity]))
            lines(p$w / p$Winf, p$g, lwd=3)
            # Plot eta_m (maybe include Wfs here?)
            abline(v=p$eta_m, lty=2, lwd=1.5)
            #  Plot selectivity axis
            axis(4, at=c(ylim.min, ylim.max + ylim.max*0.6), labels=NA, col.axis="#7E6148B2")
            mtext(c(0,100),at=c(ylim.min,  ylim.max), side=4, line=0.5) # ylim.min * 4  + ylim.max*0.6
            mtext(side=4, text="% mature individuals", line=1.2, col="#7E6148B2") # at=ylim.min * 2
            # Create logarithmic x axis
            pow <- 1:3
            ticksat <- as.vector(sapply(pow, function(p) (2:10)*10^p))
            axis(2, 10^pow, tcl=0.5, labels = NA)
            mtext(10^pow, 2, 0.5, at=10^pow)
            axis(2, ticksat, labels=NA, tcl=0.25, lwd=0, lwd.ticks=1)
            title(ylab="Growth rate (g/y)", line=1.4)
            pow <- -2:0
            ticksat <- as.vector(sapply(pow, function(p) (2:10)*10^p))
            axis(1, 10^pow, tcl=0.5, labels=NA)
            axis(1, ticksat, labels=NA, tcl=0.25, lwd=0, lwd.ticks=1)
            # Include legend with lines and relevant growth parameters
            
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
                 xlab = "",  log = "x", ylab = "",xaxt = "n", yaxt = "n", ... )
            title(xlab=expression(w/W[infinity]))
            title(ylab=expression(Mortality~(y^{-1})), line = 1.4)
            # Plot fishing mortality
            lines(p$w / p$Winf, p$psi_F * p$Fm, lwd = 3, lty = "dotted", col="#7E6148B2")
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
            # Include legend
            legend("topright", NULL, c("Natural mortality","Fishing mortality"),
                   lty=1, lwd=2, seg.len=5, bty = "n") #, "Data kernel density" ,"red" col=c("blue"),
            
            invisible(NULL)
          })

##' Returns the correlation matrix of the parameters
##'
##' 
##' @param object A \code{Parameters} object returned by \code{estimateParam}, i.e. having
##' an attribute hessian
##' @return The correlation matrix
##' @author alko
##' @rdname getCor-methods
##' @docType methods
##' @export
setGeneric("getCor", function(object) {standardGeneric("getCor")})
##' @aliases getCor,Parameters-method
##' @rdname getCor-methods
setMethod("getCor", c("Parameters"), function(object) {
  if(is.null(attr(object, "hessian"))) {
    warning("Object does not contain hessian")
    return (NULL)
  }
  cov2cor(solve(attr(object, "hessian")))
})


##' @param object a \code{\link{Parameters}} object
##'
##' @param nsim number of individuals in the simulated sample
##' @param seed the seed that is passed to \code{\link{set.seed}}.
##' @param binsize numeric, the width of the weight classes in grams
##' @param keepZeros logical, if TRUE keep bins with zero frequency. Note that the zeros after the last non-zero bin are always droped.
##'
##' @export
##' @rdname Parameters
simulate.Parameters <- function(object, nsim = 1000, seed = NULL, binsize = 5, keepZeros = TRUE, ...) {
  if(!is.null(seed)) set.seed(seed)
  p <- getParams(object)
  l <- seq(binsize / 2, p$Winf - binsize / 2, binsize)
  u <- seq(binsize + binsize / 2, p$Winf + binsize / 2, binsize)
  pr <- mapply(function(low,up){integrate(p$pdfN.approx, lower = low, upper = up)$value}, l, u)
  s <- rmultinom(1, nsim, pr)
  df <- structure(data.frame(Weight = seq(binsize, p$Winf + binsize / 2, binsize), Freq = s[,1]), binsize = binsize)
  #df <- changeBinsize2(df, binsize = binsize, keepZeros = keepZeros)
  rl <- rle(df$Freq)
  df <- head(df, nrow(df) - if(tail(rl$values, 1) == 0) tail(rl$lengths, 1) else 0)
  class(df) <-c("WeightFreq", "data.frame")
  df
}

##' @export
plot.WeightFreq <- function(x, ...) {
  plot(x, ...)
}

##' @export
hist.WeightFreq <- function(x, ..., main = "", xlab = "Weight") {
  hist(rep(x$Weight, x$Freq), xlab = xlab,  main = main, ...)
}