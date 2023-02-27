#' The negative log likelihood for given parameters and data
#'
#' \code{minimizeme} is mainly used by an optimizer (e.g. \code{optim})
#' to estimate parameters.
#'
#'
#' @param theta Numeric vector of *transformed* parameter values
#' @param data Numeric vector, or \code{data.frame} with columns named Weight and Freq. The observed values, fished individuals in grams
#' @param names String vector. Contains the names of the parameter vector theta.
#' @param fixed.names String vector. Names of constants.
#' @param fixed.vals Numeric vector. Transformed values of constants.
#' @param isSurvey Logical. If TRUE the observations are assumed to be from a survey.
#' @return Numeric scalar. The negative log likelihood for the given parameters and observations.
#' @author alko
#' @keywords optimize
#' @note If data is a list containing both sample and df, the \code{data.frame} df will be used.
#' @examples
#'
#' \dontrun{
#' ## Simulate some data with default parameter values and fishing mortality Fm = 0.3
#' sim <- simulateData3(params = s6params(Fm = 0.3))
#'
#' ## Plotting the negative log likelihood for different values of Fm
#' Fm <- seq(0.1, 1, 0.05)
#' nll <- sapply(Fm, function(x) minimizeme(theta = list(Fm = x), data = s6input(sim$df)))
#' plot(Fm, nll, type="l")
#'
#' ## Using optimise to estimate one parameter
#' wrapper <- function(x, ...) minimizeme(theta = list(Fm = x), ...)
#' est <- optimise(f = wrapper, interval=c(0.01, 1), data = s6input(sim$df))
#' est.Fm <- est$minimum
#' abline(v=est.Fm)
#' mtext(paste("Estimated Fm = ", round(est.Fm, 3)), at=est.Fm)
#' }
#'
#' @export minimizeme
#' @rdname minimizeme
minimizeme <- function(theta, data, fixed.vals = list()) {
  stopifnot(is.s6input(data))
  if (! is.null(data$surwf)) stop("Survey-data model is not implemented yet.")
  if (is.null(data$wf)) stop("No data provided")
  params <- do.call(s6params, c(theta, fixed.vals))
  return(with(getParams(params),
              sum( - data$wf$Freq * log(pdfN.approx(data$wf$Weight)) )))
}



#' Estimates parameters using the given data and maybe some parameters
#'
#' Parameter estimation minimizing the negative log likelihood, using the
#' nlminb function. \code{estimateParam} uses data from commercial catches or surveys.
#' \code{estimateMultidata} uses both sources
#'
#'
#' @param names string vector, the parameters to be estimated.
#' @param data numeric vector, or \code{data.frame} with columns Weight and Freq, or \code{list} with a numeric vector named `sample` or a \code{data.frame} named `df`.
#' Weight of individual fish (vector) or frequencies per weight class \code{data.frame}.
#' @param start numeric vector, initial values of the parameters.
#' @param lower The lower bound for the parameter estimation.
#' @param upper The upper bound for the parameter estimation.
#' @param fixed.names String vector. Names of constants.
#' @param fixed.vals Numeric vector. Transformed values of constants.
#' @param fixed.transformed Logical. If FALSE the constants are not transformed.
#' @param plotFit logical, if TRUE a plot is produced with the fited pdf and the
#' kernel density estimate of the data.
#' @param isSurvey logical, if TRUE the data are assumed to be from a survey.
#' @param verbose logical, if TRUE the estimated confidence intervals are printed.
#' @param useTMB logical, if TRUE TMB is used for parameter estimation
#' @param ... Additional named arguments passed to plotFit
#' @return \code{link{Parameters}} object, containing the estimated parameters.
#' @note If data is a list containing both sample and df, the \code{data.frame} df will be used.
#' @author alko
#' @keywords optimize
#' @examples
#'
#' ## Simulate some data
#' sam <- simulateData3(params=parameters("a", 0.5, transformed=FALSE))
#'
#' ## Estimate the a parameter and see the fitted plot
#' estimateParam(names="a", data=sam$sample, plotFit=TRUE)
#' @rdname estimateParam
#' @export estimateParam
estimateParam <-
  function(names = c("Fm", "Winf", "Wfs"),
           data=simulateData3(parameters(), samplesize=1000),
           start= rep(0.5, length(names)),
           lower=rep(-Inf, length(names)), upper=rep(Inf, length(names)),
           fixed.names=c(), fixed.vals=numeric(0), fixed.transformed = TRUE,
           plotFit=FALSE, isSurvey=FALSE, verbose=getOption("verbose"), useTMB = TRUE, ...) {
    p <- parameters()

    if(is(data, "list")) {
      if("df" %in% names(data)) {
        data <- data$df
      } else if ("sample" %in% names(data)) {
        data <- data$sample
      } else stop("`data` is a list not containing an element named sample or df.")
    }

    start[which(names == "Winf")] <-
      ifelse(is(data, "data.frame"), (max(data$Weight) + 1) / p@scaleWinf, (max(data) + 1) / p@scaleWinf)

    scales <- sapply(names, function(n) get(paste0("getscale", n))(p))

    if( ! fixed.transformed) {
      fixed.scales <- sapply(fixed.names, function(n) get(paste0("getscale", n))(p))
      fixed.vals <- log(fixed.vals / fixed.scales)
    }

    useapply <- if(require(parallel)) mclapply else lapply

    sd <- mean <- rep(0, length(names))

    estim <- nlminb(log(start), minimizeme, data=data, names=names,
                    fixed.names=fixed.names, fixed.vals=fixed.vals,
                    lower=lower, upper=upper, isSurvey=isSurvey)
    if(estim$convergence != 0) warning(estim$message)

    res <- estim$par
    h <- hessian(minimizeme, estim$par, data=data,
                 names=names, fixed.names=fixed.names,
                 fixed.vals=fixed.vals,isSurvey=isSurvey)
    s <- jacobian(minimizeme, estim$par, data=data,
                  names=names, fixed.names=fixed.names,
                  fixed.vals=fixed.vals,isSurvey=isSurvey)
    vcm <- try(solve(h))

    ci <- matrix(rep(NA, length(names)*3),ncol=3, dimnames = list(names, c("Estimate","Lower", "Upper")))
    st.er <- NA
    if( ! is(vcm, "try-error")) {
      st.er <- sqrt(diag(vcm))
      ci <- cbind(exp(res)*scales, exp(outer(1.96 * st.er, c(-1,1), '*') + res) * scales)
    }

    p <- parameters(c(names, fixed.names), c(t(simplify2array(res)), fixed.vals))

    if(plotFit) plotFit(p, data, ...)
    if(verbose) print(ci)
    return(structure(p, par=res, hessian=h, jacobian=s, st.er=st.er, ci=ci,
                     objective=estim$objective, convergence=estim$convergence,
                     nlminbMessage=estim$message, call=match.call(), version=getVersion()))
  }

##' @param surdata Same as data. Survey data.
##' @param comdata Same as data. Commercial data.
##' @rdname minimizeme
minimizemeMultidata <- function(theta, surdata, comdata, names, fixed.names=c(), fixed.vals=c())
{
  params <- parameters(c(names, fixed.names), c(theta, fixed.vals))
  return(with(getParams(params,isSurvey=TRUE), sum(-log(pdfN.approx(surdata)))) +
           with(getParams(params,isSurvey=FALSE), sum(-log(pdfN.approx(comdata)))))
}
##' @param surdata Same as data. Survey data.
##' @param comdata Same as data. Commercial data.
##' @rdname estimateParam
estimateMultidata <-
  function(names = c("Fm", "Winf", "Wfs"),
           surdata=simulateData3(parameters(), samplesize=1000, isSurvey=TRUE)$sample,
           comdata=simulateData3(parameters(), samplesize=1000)$sample,
           start= rep(0.5, length(names)), lower = rep(0.1, length(names)),
           fixed.names=c(), fixed.vals=numeric(0),
           plotFit=FALSE, ...) {
    start[which(names == "Winf")] <- (max(surdata, comdata) + 1) / parameters()@scaleWinf
    lower[which(names == "eta_F")] <- 0.0001

    useapply <- ifelse(require(parallel), mclapply, lapply)
    sd <- mean <- rep(0, length(names))


    estim <- nlminb(log(start), minimizemeMultidata, surdata=surdata, comdata=comdata,
                    names=names, fixed.names=fixed.names, fixed.vals=fixed.vals,
                    lower=log(lower))
    if(estim$convergence != 0)
      warning(estim$message)
    res <- c(estim$par) ##, estim$convergence)

    p <- parameters(c(names, fixed.names), c(t(simplify2array(res)), fixed.vals))
    if(plotFit) plotFit(p, data, ...)
    return(p)
  }

##' @export
estimate_TMB <- function(inp, n=0.75, epsilon_a=0.8, epsilon_r=0.1, A=4.47,
                         eta_m=0.25, a=0.22, Winf = NULL, sigma=NULL, u = 10,
                         sdloga = 0.7, winf.ubound = 2, Wfs = NULL,
                         verbose=options()$verbose, map=list(loga = factor(NA)),
                         random = c(), isSurvey = FALSE, eta_S = NULL, usePois = TRUE,
                         totalYield = NULL, perturbStartingVals = FALSE, ...) {
  if (is.null(df)) return(NULL)
  stopifnot(is.s6input(inp))
  if (!require(TMB)) stop("TMB is not installed! Please install and try again.")
  yrs <- names(inp$years)
  df <- df2matrix(inp$wf)
  nyrs <- ncol(df)
  DLL <- "s6modelts"

  tryer <- try({
    binsize <- attr(df,"binsize")
    if (is.null(Winf)) {
      ## Initial guess about Winf based on maximum observed
      Winf <- (nrow(df) + 2) * binsize
    } else {
      map$logWinf  <- factor(NA)
    }
    if (is.null(eta_S))  {
      ## Initial guess about eta_S
      eta_S <- which(apply(df, 1, function(x) sum(x) != 0))[[1]] * binsize / Winf
    } else {
      map$logeta_S  <- factor(NA)
    }
    if (! isSurvey) {
      map$logeta_S <- factor(NA)
    }
    if (is.null(sigma) || is.na(sigma)) {
      sigma <- if (usePois) colSums(df) else rep(0.0001, nyrs)
    } else {
      map$logSigma  <- rep(factor(NA), nyrs)
    }
    if (is.null(Wfs) || is.na(Wfs)) {
      ## Initial guess about Wfs
      Wfs <- apply(df, 2, function(x) round((which.max(x) + which(x > 0)[1]) / 2) * binsize - binsize / 2)
    } else {
      map$logWfs  <- rep(factor(NA), nyrs)
    }
    if (is.null(u)) {
      u <- 10
    } else {
      map$logu  <- factor(NA)
    }
    freq <- df
    nwc <- attr(df, "nwc")
    logFm <- rep(log(0.5), nyrs)
    data <- list(binsize = binsize, nwc = nwc, freq = freq, n = n,
                 epsilon_a = epsilon_a, epsilon_r = epsilon_r, A = A,
                 eta_m = eta_m, meanloga = log(a), sdloga = sdloga,
                 isSurvey = as.integer(isSurvey), usePois = as.integer(usePois),
                 totalYield = inp$catch$Catch)
    pars <- list(loga = log(a), logFm = logFm, logWinf = log(Winf),
                 logWfs = log(Wfs), logSigma = log(sigma),logeta_S = log(eta_S),
                 logu = log(u))
    obj <- MakeADFun(data = data, parameters = pars,  map=map, random=random, DLL = DLL)
    upper <- rep(Inf, length(obj$par))
    upper[which(names(obj$par) == "logWinf")] <- log(Winf * winf.ubound)
    obj$env$tracemgc <- verbose
    obj$env$inner.control$trace <- verbose
    obj$env$silent <- ! verbose
    if(! verbose) {
      newtonOption(obj=obj, trace=0)
      config(trace.optimize = 0, DLL=DLL)
    }
    opt <- nlminb(obj$par, obj$fn, obj$gr, upper = upper)
    sdr <- sdreport(obj, getJointPrecision = FALSE, getReportCovariance = FALSE)
    parnms <- c("Fm","Winf","Wfs", "a", "eta_S")
    vals <- sdr$value
    nms <- names(vals)
    sds <- setNames(sdr$sd, nms)
    estpars <- sapply(seq(nyrs), function(i) {
      vls <- setNames(sapply(parnms, function(x) {
        match <- vals[grepl(paste0("^", x, "$"), nms)]
        if (length(match) == 1) match else match[i]
      }), parnms)
      extrapars <- c(n, epsilon_a, epsilon_r, A, eta_m)
      names(extrapars) <- c("n", "epsilon_a", "epsilon_r", "A", "eta_m")
      do.call("s6params", as.list(c(vls, extrapars)))
    })
    Fmsy <- sapply(estpars, calcFmsy)
    SSBrel <- sapply(estpars, function(x) getParams(x)$B) /
              mapply(function(p, fmsy) getParams(s6params(Fm = fmsy, base = p))$B, estpars, Fmsy, SIMPLIFY = TRUE)
    Bexplrel <- sapply(estpars, function(x) getParams(x)$Bexpl) /
      mapply(function(p, fmsy) getParams(s6params(Fm = fmsy, base = p))$Bexpl, estpars, Fmsy, SIMPLIFY = TRUE)
    # if(nyrs == 1) {
    #   estpars <- estpars[[1]]
    #   Fmsy <- Fmsy[[1]]
    # }
    opt$convergence
 }, silent = !verbose)
  if(class(tryer) == "try-error") {
    return(tryer)
  }
  nw <- function(x) grepl(paste0("^", x, "$"), nms)
  structure(data.frame(Fm=vals[nw("Fm")], Fm_sd = sds[nw("Fm")],
                       Winf=rep(vals[nw("Winf")], nyrs), Winf_sd=rep(sds[nw("Winf")], nyrs),
                       Fmsy = Fmsy, FFmsy = vals[nw("Fm")]/Fmsy,
                       Wfs = vals[nw("Wfs")], Wfs_sd = sds[nw("Wfs")],
                       a = rep(vals[nw("a")], nyrs), eta_S = vals[nw("eta_S")],
                       sigma = vals[nw("sigma")], sigma_sd = sds[nw("sigma")],
                       u = vals[nw("u")], u_sd = sds[nw("u")],
                       R = vals[nw("R")], R_sd = sds[nw("R")],
                       Rrel = vals[nw("Rrel")], Rrel_sd = sds[nw("Rrel")],
                       Rp = vals[nw("Rp")], Rp_sd = sds[nw("Rp")],
                       rmax = vals[nw("rmax")], rmax_sd = sds[nw("rmax")],
                       Y = vals[nw("Y")], Y_sd = sds[nw("Y")],
                       ssb = vals[nw("ssb")], ssb_sd = sds[nw("ssb")],
                       Bexpl = vals[nw("Bexpl")], Bexpl_sd = sds[nw("Bexpl")],
                       ssbrel = SSBrel, Bexplrel = Bexplrel,
                       row.names=yrs),
            obj=obj, opt=opt, sdr = sdr, estpars=estpars, inp = inp)
}

