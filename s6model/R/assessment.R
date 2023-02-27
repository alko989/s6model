##' Class s6modelResults
##'
##' @title s6modelResults class
##' @exportClass s6modelResults
##'@name s6modelResults
#' @aliases s6modelResults-class
#' @rdname s6modelResults
#' @export
setClass("s6modelResults")

#' Fit weight-length relationship
#'
#' @param df \code{data.frame} or \code{DATRASraw} object containing length and weight data
#' @param colname.weight chr, the name of the column containing weight data. Only used if df is \code{data.frame}
#' @param colname.length chr, the name of the column containing length data. Only used if df is \code{data.frame}
#' @param plotFit logical, if TRUE plot the fitted line along with the input data
#' @param mindata integer, the minumum data points. If less data are available, default paramter values are returned with a warning.
#' @param ... additional arguments passed to plot
#'
#' @return A list with a and b parameters and data points used (n)
#' @export
#'
#' @examples
#' ## Simulated data
#' a <- 0.00058
#' b <- 3.123
#' df <- data.frame(lengths = 50:100,  weights =  a * c(50:100) ^ b + rnorm(51))
#' 
#' ## Fit weight-length relationship and compare with simulated values of a and b
#' s6model:::fitWL(df, colname.weight = "we", colname.length = "le")
#' ## a: 0.00057 (sim: 0.00058)
#' ## b: 3.128   (sim: 3.123)
#' 
#' @note \code{fitWL} return an error if data from more than one species are contained in the DATRASraw object.
#' @note If mindata is 4 or lower the minimum acceptable amount of data points is set to 4.
fitWL <- function(df, colname.weight = "Weight", colname.length = "Length", plotFit = FALSE, mindata = 0, ...) {
  if(is(df, "data.frame")) {
    mindata <- max(4,mindata)
    nms <- names(df)
    wlcols <- c(pmatch(colname.weight, nms), pmatch(colname.length, nms))
    df <- setNames(df[wlcols], c("Weight", "Length"))
    w <- which(df$Weight > 0 & df$Length > 0, ! is.na(df$Weight) & ! is.na(df$Length))
    df <- df[w, ]
    if(mindata > 0 & nrow(df) < mindata) {
      if(plotFit) {
        plot(1, axes = FALSE, type = "n", xlab = "", ylab = "", ...)
        text(1,1, adj = 0.5, "No weight-length data")
      }
      warning("Not enough data. Default parameter values are returned")
      return(list(a = 0.01, b = 3, n = nrow(df)))
    }
    fit <- lm(log(Weight) ~ log(Length), data = df) 
    a <- as.numeric(exp(fit$coefficients[1]))
    b <- as.numeric(fit$coefficients[2])
    res <- list(a = a, b = b, n = nrow(df))
  } else if(is(df, "DATRASraw")) {
    if(nlevels(df[[3]]$Species) != 1) {
      stop("The DATRASraw object contains data for more than one species")
    }
    if(levels(df[[3]]$Species) != levels(df[[1]]$Species)) {
      warning("The CA and HL parts of the DATRASraw object do not have the same species")
    }
    return(fitWL(df[[1]], colname.weight = "IndWgt", colname.length = "LngtCm", plotFit = plotFit))
  }
  if(plotFit) {
    xlim <- range(0, df$Length)
    plot(Weight ~ Length, data = df, pch = ".", cex = 1.1, xlim = xlim, ...)
    newdat <- data.frame(Length = seq(0, max(df$Length), length.out = 100))
    lines(newdat$Length, exp(predict(fit, newdata = newdat)), col =2, lwd = 2)
    ss <- summary(fit)
    mtext(side = 3, adj = 0, at = 0, bquote(W == .(round(res$a,4))* L^.(round(res$b,2))), line = -2)
    mtext(side = 3, adj = 0, at = 0, bquote(paste(adj.R^2, "=", .(round(ss$adj.r.squared, 2)))), line = -3.2)
    mtext(side = 3, adj = 0, at = 0, bquote(paste(n, "=", .(nrow(df)))), line = -4.4)
  }
  res
}

l2w <- function(l, a, b) {
  a * l ^ b
}

addWeight <- function(df, a, b, lengthcol = "Length") {
  if(is.null(df$Weight))
    df$Weight <- l2w(l = df[[lengthcol]], a = a, b = b)
  else
    warning("Column `Weight` exists, the original data frame `df` is returned")
  df
}


getalim <- function (p) {
  if(!is.s6params(p)) { 
    stop("p should be a s6params object") ## return(0.8)
  }
  optimize(function(x) {
    p$a <- x
    getParams(p, optim.Rrel = TRUE, FF = 0)^2
  }, c(0,2) )$minimum
}


getCI <- function (inputData, ests, a.mean, a.sd, same.as = TRUE, nsample, winf.ubound, yield, probs, Winf, u, ...) {
  if(is(ests, "Parameters")) {
    ests <- list(ests)
  }
  aplfun <- if(require(parallel)) mclapply else lapply
  maplfun <- if(require(parallel)) mcmapply else mapply
  r <- function(x) round(x, 2)
  if(same.as) {
    alim <- getalim(meanParameters(ests))  
    as <- rtrunc(nsample, spec ="lnorm",  meanlog = log(a.mean), sdlog = a.sd, 
                 b = alim)
  }
  ci <- lapply(seq(along.with = inputData), function(i) {
    if(!same.as) {
      alim <- getalim(ests[[i]])
      as <- rtrunc(nsample, spec ="lnorm", meanlog = log(a.mean), sdlog = a.sd, 
                   a = -Inf, b = alim)
    }
    Winf <- if(is.null(ests[[i]])) NULL else getWinf(ests[[i]])
    results <- aplfun(as, function(a) {
      estimate_TMB(inputData[[i]], a = a, Winf = Winf, totalYield = yield[i],
                   u = u, ...)
    })
    reps <- results
    asused <- as.list(as)
    err <- sapply(results, function(x) is(x, "try-error"))
    errmsg <- reps[err]
    reps[err] <- NULL
    asused[err] <- NULL
    notConv <- sapply(reps, function(x) attr(x, "opt")$convergence == 1)
    reps[notConv] <- NULL
    asused[notConv] <- NULL
    verySmall <- sapply(reps, function(x) x$Fm < 1e-5)
    reps[verySmall] <- NULL
    asused[verySmall] <- NULL
    negRrel <- sapply(reps, function(x) x$Rrel < 0)
    reps[negRrel] <- NULL
    asused[negRrel] <- NULL
    repsdf <- do.call(rbind.data.frame, reps)
    nrep <- nrow(repsdf)
    n <- function(x) if(length(x) > 0) sum(x) else 0
    resdf <- as.data.frame(apply(repsdf, 2, quantile, probs=probs, na.rm=TRUE))
    structure(resdf, results = results, alim=alim, as = as, nrep = nrep,
              notConv = n(notConv) , nVerySmall = n(verySmall), nerr = n(err), 
              err = err, errmsg = errmsg, pointests = repsdf, asused = asused,
              nnegRrel = n(negRrel))
  })
  alims <- sapply(ci, function(x) attr(x, "alim"))
  reps <- sapply(ci, function(x) attr(x, "nrep"))
  nerr <- sapply(ci, function(x) attr(x, "nerr"))
  nnegRrel <- sapply(ci, function(x) attr(x, "nnegRrel"))
  nVerySmall <- sapply(ci, function(x) attr(x, "nVerySmall"))
  err <- sapply(ci, function(x) attr(x, "err"))
  errmsg <- sapply(ci, function(x) attr(x, "errmsg"))
  notConv <- sapply(ci, function(x) attr(x, "notConv"))
  pointests <-  sapply(ci, function(x) attr(x, "pointests"))
  ## allResults <- lapply(ci, function(x) attr(x, "results"))
  as <-  lapply(ci, function(x) attr(x, "as"))
  asused <-  lapply(ci, function(x) attr(x, "asused"))
  nms <- names(ci[[1]])
  ci <- lapply(nms, function(nm) lapply(ci, function(d) d[[nm]]))
  ci <- lapply(ci, function(yy) {
    for(w in which(sapply(yy, is.null))) {
      yy[[w]] <- rep(NA, length(probs))
    }
    setNames(do.call(cbind.data.frame, yy), names(inputData))
  })
  ci <- setNames(ci, nms)
  structure(ci, alims = alims, as = unlist(as), asused = unlist(asused),
            nVerySmall = nVerySmall, notConv = notConv, nnegRrel = nnegRrel,
            pointests = pointests, err = nerr, err = err, errmsg = errmsg, 
            reps = reps)
}


df2matrix <- function(df){
  maxrow <- max(sapply(df, nrow))
  structure(sapply(df, function(x) c(x$Freq, rep(0, maxrow - nrow(x))) ),
            nwc = maxrow, binsize = attr(df[[1]], "binsize"))  
}

#' Makes assessment of a fish stock given weight frequency data for several years
#'
#' @param inputData list of data.frames, each data.frame has columns \code{Weight} and  \code{Freq} and attribute \code{binsize}
#' @param yield numeric, the total yearly catch or landings in kg. Use NULL if not known.
#' @param a.mean numeric, physiological mortality.
#' @param a.sd numeric, the standard deviation (log domain) of the log-normal distribution of physiological mortality.
#' @param nsample integer, number of repetitions for uncertainty etimation. If zero no uncertainty is estimated.
#' @param same.as logical, if TRUE use the same random values of physiological mortality for each year.
#' @param seed numeric, the random number generator seed.
#' @param u numeric, the selectivity steepness parameter.
#' @param sigma numeric, if NULL the parameter is estimated, otherwise a constant is used, see Details.
#' @param binsize numeric, span of weight classes in grams.
#' @param winf.ubound numeric, the upper bound of asymptotic weight. It is a multiplier of the maximum observed weight.
#' @param equalWinf logical, if TRUE estimate one asymptotic weight for all years, if FALSE estimate one for each year.
#' @param probs numeric vector of probabilites with values in [0,1] for the uncertainty sample quantiles.
#' @param dirout Output directory
#' @param fnout Output file name
#' @param ... Arguments passed to \code{estimate_TMB} that does the estimation
#'
#' @return object of class s6modelResults which is a data.frame with parameter estimates with attributes \describe{
#' \item{CI}{confidence levels as sample quantiles}
#' \item{obj}{the TMB object of the default value estimation}
#' \item{opt}{the optimization output from \code{nlminb}}
#' \item{opts}{the input options used in the run}
#' \item{Results}{all results using default parameters}
#' \item{seed}{the random number gernerator used}
#' \item{timeToCompletion}{difftime object with the time needed for the estimation}
#' \item{version}{the \code{s6model}version that produced the results}
#' }  
#' 
#' @details \itemize{
#' \item{sigma paramter} for more information see the cited paper.
#' }
#' 
#' @export
makeAssessment <- 
  function(inputData, yield = NULL, a.mean = 0.22, a.sd = 0.7, nsample = 100, 
           same.as = TRUE, u = 10, sigma = NULL, binsize = NULL, 
           winf.ubound = 2, equalWinf = TRUE, probs = seq(0, 1, 0.01), 
           seed = as.integer(rnorm(1, 1000, 100)), dirout = "results",
           fnout = format(Sys.time(),"results_%Y%m%d_%H%M.RData"), ...) {
  set.seed(seed)
  haveYield <- TRUE
  if(is.null(yield)) {
    yield <- rep(1, length(inputData))
    haveYield <- FALSE
  }
  sigma <- if(is.null(sigma) || is.na(sigma)) {
    rep(NA, length(inputData)) 
  } else {
    sapply(inputData, function(x) mean(rle(x$Freq)$lengths) * sum(x$Freq))
  }
  starting <- Sys.time()
  if(! is.null(binsize)) {
    inputData <- changeBinsize2(inputData, binsize = binsize)
  }
  if(equalWinf) {
    ests <- try(
      estimate_TMB(inputData, DLL = "s6modelts", totalYield = yield, u = u,
                   sigma = sigma, a = a.mean, winf.ubound = winf.ubound, ...))
    estpars <- attr(ests, "estpars")
    res <- ests
  } else {
    ests <- mapply(function(x, y, s) estimate_TMB(x, a=a.mean, winf.ubound = winf.ubound,
                                                  totalYield = y, sigma = s, u = u, DLL = "s6model", ...),
                   inputData, yield, sigma, SIMPLIFY = FALSE)
  res <- lapply(ests, function(x) if(class(x)== "try-error") rep(NA, 30) else x[1:30] )
    res <- do.call(rbind.data.frame, res)
    row.names(res) <- names(inputData)
    estpars <- lapply(ests, function(x) attr(x, "estpars"))
    attr(res, "estpars") <- estpars
  }
  if(is(res, "try-error")) {
    cat("There was an error whith the estimation")
    return(res)
  }  
  if(a.sd > 0 & nsample > 1) {
    attr(res, "CI") <- getCI(inputData = inputData, ests = estpars, a.mean = a.mean, a.sd = a.sd, u = u,
                             same.as = same.as, nsample = nsample, winf.ubound = winf.ubound, yield = yield, probs = probs, ...)
  } 
  opts <- list(...)
  res <- structure(res, Results = ests, version = getVersion(), timeToCompletion = Sys.time() - starting, haveYield = haveYield,
                   seed = seed, opts = list(tmbopts = c(opts, a.mean = a.mean, a.sd = a.sd, 
                                                        winf.ubound = winf.ubound, yield = yield, u = u), probs = probs))
  if( ! is(res,"try-error")) { 
    class(res) <- c("s6modelResults", class(res))
  }
  if(!is.null(fnout)) {
    dir.create(dirout, showWarnings = FALSE)
    out <- file.path(dirout, fnout)
    save(res, file = out)
    cat("Results are saved in:", out, "\n")
  }
  res
}

##' @export
##' @rdname s6modelResults
print.s6modelResults <- function(x, ...) {
  n <- dim(x)
  if(dim(x)[1] == 0) {
    cat("Object with no results")
  } else {
    print.data.frame(x, ...)
  }
  cat("\nResults produced by:", formatVersion(attr(x, "version")), "\n")
  invisible(x)
}

#' @rdname s6modelResults
#' @export
residuals.s6modelResults <- function(object, ...) {
  rprt <- attr(object, "obj")$env$report()
  w <- which(rprt$freq > 0)
  rprt$residuals[w]
}


makeShading <- function(x, ylow, yhigh, col = grey(0.8), alpha = 1, ...) {
  xs <- c(x, rev(x))
  ys <- c(ylow, rev(yhigh))
  notna <- ! is.na(ys)
  col <- col2rgb(col) / 255
  col <- rgb(col[1], col[2], col[3], alpha)
  polygon(xs[notna], ys[notna], col = col, border = NA)
}

addConfidenceShading <- 
  function(x, y, ..., probs = c(0.025, 0.975), exclude = 0, 
           what = "FFmsy", grey.intensity = 1.5, col.ci = 1, alpha = 1,
           addMedian = FALSE, col.median = "white", lty.median = 2, lwd.median = 2) {
    if(is(y, "s6modelResults")) {
      r <- attr(y, "Results")
      w <- sapply(r, function(rr) {
        if(is(rr, "try-error")) return(rep(NA, 2))
        sdr <- attr(rr, "sdr")
        id <- which(names(sdr$value) == ifelse(what == "FFmsy", "Fm", what))
        c(sdr$value[id] + qnorm(probs) * sdr$sd[id])
      })
      if(what == "FFmsy") {
        w <- sweep(w, 2, y$Fm / y$FFmsy, "/")
      }
      makeShading(x, w[1,], w[2, ],  col = col.ci, alpha = alpha)
    } else if (is(y, "data.frame")){
      d <- nrow(y) - 1
      for(i in seq(1 + exclude, d / 2)) {
        gr <- max(min(1, (1 - i / ((d - (2 * exclude))) / grey.intensity )), 0)
        makeShading(x, y[i, ], y[d - i, ], col=grey(gr), alpha = alpha)
      }
      if(addMedian) {
        lines(x, y[d/2 + 1, ], col = col.median, lty = lty.median, lwd=lwd.median)
      }
    } else {
      stop("Only data.frame or s6modelResults are accepted as input to addConfidenceShading")
    }
  }

getYears <- function(x) {
  as.numeric(regmatches(x, regexpr("[0-9]+", x)))
}

constrFilename <- function(stock, a, a.sd, winf.ubound, Winf, sigma, usePois, u,
                           aggryrs, nsample, includeUncertainty = nsample > 1, equalWinf, binsize, ...) {
  estimateWinf <- is.null(Winf)
  paste0(stock, "_a=", a, 
         "_asd=", a.sd, 
         if(estimateWinf) paste0("_estWinf_winfUbound=", winf.ubound) else paste0("_fixWinf=", Winf), 
         if(is.null(sigma)) "_estSigma" else paste0("_sigma=", sigma), 
         if(usePois) "_usePoison" else "_useGauss", "_aggryrs=", aggryrs,  
         if(is.null(u)) "_estu" else "_fixu=10",
         "_nsample=", if(includeUncertainty) nsample else 1, 
         if(equalWinf) "_equalWinf" else "_difWinf", "_binsize=", binsize, "_", format(Sys.time(), format = "%Y%m%d_%H%M%S"),
         ".RData")
}

##' @export
##' @rdname s6modelResults
plot.s6modelResults <- function(x, ..., what = "FFmsy", use.rownames = TRUE, 
                                years = NULL, xlab = NULL, ylab = NULL, 
                                ylim = NULL, addDefault = FALSE, col.def = "white",
                                addhline = 1, col.hline = 1, lty.hline = 2,
                                addMedian = TRUE, col.median = "white", lty.median = 2, lwd.median = 2,
                                cex.ver = 0.7, version = FALSE, xaxs = "i", yaxs = "i",
                                ci = c("bootstrap", "estimated", "none"), grey.intensity = 1.5,
                                exclude = 0, mult = 1, add = FALSE, alpha = 1) {
  ci <- match.arg(ci)
  ys <- x[[what]] / mult
  CI <- attr(x, "CI")
  if(ci == "bootstrap" & ! is.null(CI)) {
    cidf <- CI[[what]] / mult
    if(is.null(cidf)) {
      warning("There is no confidence interval for ", what, " in the results. Only the default line will be plotted")
      ci <- "none"
      addDefault <- TRUE
    } else {
      d <- nrow(cidf) - 1
      median <- cidf[d/2 + 1, ]
    }
  }
  yl <- switch(what, FFmsy = expression(F/F[msy]), Fm = expression(F~(y^-1)), 
               Winf = expression(W[infinity]~(g)), 
               Wfs = "50% retainment size (g)", 
               ssb = {
                 haveYield <- attr(x, "haveYield")
                 if(is.null(haveYield) || haveYield) {
                   paste0("SSB (", ifelse(mult > 1, paste0(mult, " "), "") ,"t)")
                 } else { 
                   cidf <- data.frame(mapply(`/`, cidf, x$rmax))
                   bquote(B[SS]/R[max]~(g*y))
                 }
               },
               ssbrel = expression(B[SS]/B[SS]^{msy}),
               R = "Recruitment",
               stop("Unidentified `what` argument. Please select one of FFmsy, Fm, Winf, ssb, ssbrel, Wfs, or R"))
  ylab <- if(is.null(ylab)) yl else ylab
  xlab <- if(is.null(xlab)) "Year" else xlab
  xs <- seq(dim(x)[1])
  if(use.rownames) {
    xs <- getYears(rownames(x))
  }
  if(! is.null(years) && ! use.rownames) {
    xs <- years
  }
  if(! add) {
    ## Calculate y-axis range if ylim is not provided
    if(is.null(ylim)) {
      if(is.null(list(...)$log) || (! grepl("y", list(...)$log))) ylim <- c(0,0)
      if(addMedian) ylim <- range(ylim, median, na.rm = TRUE)
      if(addDefault) ylim <- range(ylim, ys, na.rm = TRUE)
      if(ci == "bootstrap") ylim <- range(ylim, cidf[-c(1:exclude, (nrow(cidf) - 5 + 1):nrow(cidf)), ], na.rm = TRUE)
    }    
    plot(xs, ys, type="n", ylim=ylim, xlab = "", ylab = "", xaxs = xaxs, yaxs = yaxs, ...)  
  }
  
  title(xlab = xlab, ylab=ylab, line=2)
  if(ci == "bootstrap") {
    if( ! is.null(cidf)) {
      addConfidenceShading(xs, cidf, grey.intensity = grey.intensity, 
                           addMedian = addMedian, col.median = col.median, 
                           lty.median = lty.median, lwd.median = lwd.median,
                           exclude = exclude, alpha = alpha)
    }
  } else if(ci == "estimated") {
    addConfidenceShading(xs, x, what = what)  
  } else if(ci == "none") {
    ## Nothing to see here, move along  
  } else stop(ci, " is not recognized confidence interval")
  if(addMedian) {
    lines(xs, median, col = col.median, lty = lty.median, lwd=lwd.median)
  }  
  if(addDefault){
    lines(xs, ys, lty=2, lwd=3, col= col.def)
  }
  if( ! is.null(addhline)) {
    abline(h=addhline, lwd=2, col=col.hline, lty=lty.hline)
  }
  if(version) {
    addVersion(attr(x, "version"), cex = cex.ver, lengthSHA = 6)
  }
  box()
}
