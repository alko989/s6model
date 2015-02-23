##' @title s6modelResults class
##' @exportClass s6modelResults
##'@name s6modelResults
#' @aliases s6modelResults-class
#' @rdname s6modelResults
#' @export
setClass("s6modelResults")

fitWL <- function(df, colname.weight = "Weight", colname.length = "Length") {
  if(is(df, "data.frame")) {
    df <- setNames(df[c(colname.weight, colname.length)], c("Weight", "Length"))
    w <- which(df$Weight > 0 & df$Length > 0, ! is.na(df$Weight) & ! is.na(df$Length))
    fit <- lm(log(Weight) ~ log(Length), data = df[w, ]) 
    a <- as.numeric(exp(fit$coefficients[1]))
    b <- as.numeric(fit$coefficients[2])
    res <- list(a = a, b = b)
  } else if(is(df, "DATRASraw")) {
    if(nlevels(df[[3]]$Species) != 1) {
      stop("The DATRASraw object contains data for more than one species")
    }
    if(levels(df[[3]]$Species) != levels(df[[1]]$Species)) {
      warning("The CA and HL parts of the DATRASraw object do not have the same species")
    }
    res <- fitWL(df[[1]], colname.weight = "IndWgt", colname.length = "LngtCm")
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
  if(!is(p, "Parameters")) { 
    return(0.8)
  }
  optimize(function(x) {
    getParams(parameters("a", x, FALSE, base = p), optim.Rrel = TRUE, FF = 0)^2
  }, c(0,2) )$minimum
}


getCI <- function (inputData, ests, a.mean, a.sd, nsample, winf.ubound, yield, probs, ...) {
  aplfun <- if(require(parallel)) mclapply else lapply
  maplfun <- if(require(parallel)) mcmapply else mapply
  r <- function(x) round(x, 2)
  ci <- lapply(seq(along.with = inputData), function(i) {
    alim <- getalim(ests[[i]])
    ## cat(paste0("Note: physiological mortality a ~ truncLogNorm(", 
    ##           r(log(a.mean)), ", ", r(a.sd), ", ubound = ", r(alim), ")\n"))
    as <- rtrunc(nsample, spec ="lnorm",  meanlog = log(a.mean), sdlog = a.sd, b = alim)
    results <- aplfun(as, function(a) estimate_TMB(inputData[[i]], a = a, winf.ubound = winf.ubound, 
                                                totalYield = yield[i], ...))
    reps <- results
    err <- sapply(results, function(x) is(x, "try-error"))
    errmsg <- reps[err]
    reps[err] <- NULL
    nerr <- sum(err)
    notConv <- sapply(reps, function(x) attr(x, "opt")$convergence == 1)
    reps[notConv] <- NULL
    repsdf <- do.call(rbind.data.frame, reps)
    nrep <- nrow(repsdf)
    structure(as.data.frame(apply(repsdf, 2, quantile, probs = probs, na.rm = TRUE)), results = results,
              alim=alim, as = as, nrep = nrep, notConv = if(length(notConv) > 0) sum(notConv) else 0, nerr = nerr, err = err, errmsg = errmsg)
  })
  alims <- sapply(ci, function(x) attr(x, "alim"))
  reps <- sapply(ci, function(x) attr(x, "nrep"))
  nerr <- sapply(ci, function(x) attr(x, "nerr"))
  err <- sapply(ci, function(x) attr(x, "err"))
  errmsg <- sapply(ci, function(x) attr(x, "errmsg"))
  notConv <- sapply(ci, function(x) attr(x, "notConv"))
  ## allResults <- lapply(ci, function(x) attr(x, "results"))
  as <-  lapply(ci, function(x) attr(x, "as"))
  nms <- names(ci[[1]])
  ci <- lapply(nms, function(nm) lapply(ci, function(d) d[[nm]]))
  ci <- lapply(ci, function(yy) {
    for(w in which(sapply(yy, is.null))) {
      yy[[w]] <- rep(NA, length(probs))
    }
    setNames(do.call(cbind.data.frame, yy), names(inputData))
  })
  ci <- setNames(ci, nms)
  structure(ci, alims = alims, as = as, reps = reps, notConv = notConv, nerr = nerr, 
            err = err, errmsg = errmsg) ##, allResults = allResults)
}

##' @export
makeAssessment <- function(inputData, a.mean = 0.27, a.sd = 0.89, nsample = 100, 
                           probs = seq(0, 1, 0.01), winf.ubound = 2, 
                           dirout = "results", yield = NULL, seed = as.integer(rnorm(1, 1000, 100)),
                           fnout = format.Date(Sys.time(), "results_%Y%m%d_%H%M.RData"), sigma = NULL, ...) {
  set.seed(seed)
  timeToCompletion <- system.time({
  if(is.null(yield)) yield <- rep(0.0001, length(inputData))
  sigma <- if(is.null(sigma) || is.na(sigma)) rep(NA, length(inputData)) else sapply(inputData, function(x) mean(rle(x$Freq)$lengths) * sum(x$Freq))
  ests <- mapply(function(x, y, s) estimate_TMB(x, a=a.mean, winf.ubound = winf.ubound, totalYield = y, sigma = s, ...),
                 inputData, yield, sigma, SIMPLIFY = FALSE)
  res <- lapply(ests, function(x) if(class(x)== "try-error") rep(NA, 15) else x[1:15] )
  res <- do.call(rbind.data.frame, res)
  row.names(res) <- names(inputData)
  if(a.sd > 0 & nsample > 1) {
    attr(res, "CI") <- getCI(inputData, ests, a.mean, a.sd, nsample, winf.ubound, yield, probs, ...)
  } 
  })
  opts <- list(...)
  res <- structure(res, Results = ests, version = getVersion(), timeToCompletion = timeToCompletion, seed = seed,
                   opts = list(tmbopts = c(opts, a.mean = a.mean, a.sd = a.sd, 
                                           winf.ubound = winf.ubound, yield = yield), probs = probs))
  class(res) <- c("s6modelResults", class(res))
  if(!is.null(fnout)) {
    dir.create(dirout, showWarnings = FALSE)
    out <- file.path(dirout, fnout)
    save(res, file = out)
    cat("Results are saved in:", out)
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

makeShading <- function(x, ylow, yhigh, col = grey(0.8), ...) {
  xs <- c(x, rev(x))
  ys <- c(ylow, rev(yhigh))
  notna <- ! is.na(ys)
  polygon(xs[notna], ys[notna], col = col, border = NA)
}

addConfidenceShading <- 
  function(x, y, ..., probs = c(0.025, 0.975), exclude = 0, 
           what = "FFmsy", grey.intensity = 1.5, 
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
      makeShading(x, w[1,], w[2, ],  col = "lightgrey")
    } else if (is(y, "data.frame")){
      d <- nrow(y) - 1
      for(i in seq(1 + exclude, d / 2)) {
        gr <- max(min(1, (1 - i / ((d - (2 * exclude))) / grey.intensity )), 0)
        makeShading(x, y[i, ], y[d - i, ], col=grey(gr) )
      }
      if(addMedian) {
        lines(x, y[d/2 + 1, ], col = col.median, lty = lty.median, lwd=lwd.median)
      }
    } else {
      stop("Only data.frame or s6modelResults are accepted as input by addConfidenceShading")
    }
  }

##' @export
##' @rdname s6modelResults
plot.s6modelResults <- function(x, ..., what = "FFmsy", use.rownames = TRUE, 
                                years = NULL, xlab = NULL, ylab = NULL, 
                                ylim = NULL, addDefault = FALSE, col.def = "white",
                                addhline = 1, col.hline = 1, lty.hline = 2,
                                addMedian = FALSE, col.median = "white", lty.median = 2, lwd.median = 2,
                                cex.ver = 0.7, version = FALSE, xaxs = "i", yaxs = "i",
                                ci = c("bootstrap", "estimated", "none"), grey.intensity = 1.5,
                                exclude = 0, mult = 1) {
  yl <- switch(what, FFmsy = expression(F/F[msy]), Fm = expression(F~(y^-1)), 
               Winf = expression(W[infinity]~(g)), 
               Wfs = "50% retainment size (g)", 
               ssb = paste0("SSB (", ifelse(mult > 1, paste0(mult, " "), "") ,"t)"),
               stop("Unidentified `what` argument. Please select one of FFmsy, Fm, Winf, ssb, or Wfs"))
  ylab <- if(is.null(ylab)) yl else ylab
  xlab <- if(is.null(xlab)) "Year" else xlab
  xs <- seq(dim(x)[1])
  if(use.rownames) {
    xs <- as.numeric(regmatches(row.names(x), regexpr("[0-9]+", row.names(x))))
  }
  if(! is.null(years) & ! use.rownames) {
    xs <- years
  }
  ys <- x[[what]] / mult
  cidf <- attr(x, "CI")[[what]] / mult
  ylim <- if(is.null(ylim)) range(ys, cidf, na.rm = TRUE) else ylim
  
  plot(xs, ys, type="n", ylim=ylim, xlab = "", ylab = "", xaxs = xaxs, yaxs = yaxs, ...)
  title(xlab = xlab, ylab=ylab, line=2)
  if("bootstrap" %in% ci) {
    if( ! is.null(cidf)) {
      addConfidenceShading(xs, cidf, grey.intensity = grey.intensity, 
                           addMedian = addMedian, col.median = col.median, 
                           lty.median = lty.median, lwd.median = lwd.median,
                           exclude = exclude)
    }
  } else if("estimated" %in% ci) {
    addConfidenceShading(xs, x, what = what)  
  } else if("none" %in% ci) {
    col.def <- 1
    addDefault <- TRUE
  } else stop(ci, " is not recognized confidence interval")
  
  if(addDefault){
    lines(xs, ys, lty=2, lwd=3, col= col.def)
  }
  if(!is.null(addhline)) {
    abline(h=addhline, lwd=2, col=col.hline, lty=lty.hline)
  }
  if(version) {
    addVersion(attr(x, "version"), cex = cex.ver, lengthSHA = 6)
  }
  box()
}

##' @export 
addIces<- function(stock, icesfile = "~/Work/mainCode/R/SecondPaper/ICES/ICES-cod.RData",
                   col="darkgrey", lwd=2, lty = c(2,1,2), what = "FFmsy", mult = 1, ...) {
  load(icesfile)
  ices <- ices.cod[[stock]]
  nms <- tolower(names(ices))
  if(what %in% c("FFmsy", "ffmsy"))
    matplot(ices$Year, ices[ , c("high_F", "F","low_F")] / fmsy(ices) / mult, 
            add=TRUE, col=col, lwd = lwd, lty = lty, type="l")
  else {
    n <- pmatch(what, nms)
    h <- pmatch(paste0("high_", what), nms)
    l <- pmatch(paste0("low_", what), nms)
    matplot(ices$Year, ices[ , c(h, n, l)] / mult, 
            add=TRUE, col=col, lwd = lwd, lty = lty, type="l")
  }
}