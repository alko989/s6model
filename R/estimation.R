
#' Estimate fishing mortality and stock status
#' 
#' @description Function that takes weight frequency distributions and estimates 
#' the stock status (F/Fmsy), the asymptotic weight and the 50\% retainmen size.
#'
#' @param inp \code{\link{s6input}} object containing the input data
#' @param params \code{\link{s6params}} object containing the model parameters
#' @param winf.ubound numeric, upper bound for Winf
#' @param est.Winf logical, the est.* arguments are used to estimate (TRUE) or fix (FALSE) model parameters.
#' @param est.Wfs logical
#' @param est.etaS logical
#' @param est.a logical
#' @param est.sigma logical
#' @param est.u logical
#' @param verbose logical, if \code{TRUE} print TMB output
#' @param random character vector, select which parameters should be random effects
#' @param make.guess logical, if \code{TRUE} a guess is used to initialise parameter values, otherwise the values from params is used.
#' @param usePois logical, if \code{TRUE} (default) the Poisson estimation is used, otherwise a Gaussian estimation is used.
#' @param sdloga numeric, the standard deviation of the prior of physiological mortality (a).
#' @param ... named arguments passed to s6params
#'
#' @details The upper bound of Winf is set using \code{winf.ubound}. The 
#' bound then is equal to the maximum observed weight in the data times 
#' \code{winf.ubound}.
#'
#' @return A data.frame with all parameter estimates and derived quantities and their uncertainty.
#' @export
#'
estimate <- function(inp, params = s6params(...),
                     winf.ubound = 2,
                     est.Winf = TRUE, est.Wfs = TRUE, est.etaS = inp@isSurvey, 
                     est.a = FALSE, est.sigma = TRUE, est.u = FALSE,
                     verbose = FALSE, random = c(), make.guess = TRUE,
                     usePois = TRUE, sdloga = 0.7, ...) {
  if (is.null(inp)) return(NULL)
  if ( ! is(inp, "s6input")) stop(sprintf("Function `estimate` expects a `s6input` object, a %s was provided.", is(inp)))
  dat <- df2matrix(if(inp@isSurvey) inp@surWF else inp@wf)
  binsize <- attr(dat,"binsize")
  nyrs <- ncol(dat)
  map <- list()
  p <- as.list(params)
  tryer <- try({
    if (! est.a) {
      map$loga <- factor(NA)
    }
    if (est.Winf) {
      if (make.guess) {
        p$Winf <- (nrow(dat) + 2) * binsize
      }
    } else {
      map$logWinf  <- factor(NA)
    }
    if (est.etaS & inp@isSurvey)  {
      if (make.guess) {
        p$eta_S <- which(apply(dat, 1, function(x) sum(x) != 0))[[1]] * binsize / p$Winf
      }
    } else {
      map$logeta_S  <- factor(NA)
    }
    if (est.sigma) {
      sigma <- if(usePois) colSums(dat) else rep(0.0001, nyrs)
    } else {
      map$logSigma  <- rep(factor(NA), nyrs)
    }
    if (est.Wfs)  {
      if (make.guess) {
        p$Wfs <- apply(dat, 2, function(x) round((which.max(x) + which(x > 0)[1]) / 2) * binsize - binsize / 2)
      }
    } else {
      map$logWfs  <- rep(factor(NA), nyrs)
    }
    if (est.u)  {
      if (make.guess) {
        p$u <- 10
      }
    } else {
      map$logu  <- factor(NA)
    }
    nwc <- attr(dat, "nwc")
    logFm <- rep(log(0.5), nyrs)
    data <- list(binsize=binsize, nwc=nwc, freq=dat, n=p$n, epsilon_a=p$epsilon_a,
                 epsilon_r=p$epsilon_r, A=p$A, eta_m=p$eta_m, meanloga = log(p$a), 
                 sdloga = sdloga, isSurvey = as.integer(inp@isSurvey),
                 usePois = as.integer(usePois), Catch = inp@catch)
    pars <- list(loga = log(p$a), logFm = logFm, logWinf = log(p$Winf),
                 logWfs = log(p$Wfs), logSigma = log(sigma), logeta_S = log(p$eta_S), 
                 logu = log(p$u))
    obj <- MakeADFun(data = data, parameters = pars, map=map, 
                     random=random, DLL = "s6model")
    upper <- rep(Inf, length(obj$par))
    upper[names(obj$par) == "logWinf"] <- log(p$Winf * winf.ubound)
    obj$env$tracemgc <- verbose
    obj$env$inner.control$trace <- verbose
    obj$env$silent <- ! verbose
    if(! verbose) {
      newtonOption(obj=obj, trace=0)
      config(trace.optimize = 0, DLL="s6model")
    }
    opt <- nlminb(obj$par, obj$fn, obj$gr, upper = upper)
    sdr <- sdreport(obj)
    parnms <- c("Fm","Winf","Wfs", "a", "eta_S")
    vals <- sdr$value
    nms <- names(vals)
    sds <- setNames(sdr$sd, nms)
    estpars <- sapply(seq(nyrs), function(i) {
      vls <- sapply(parnms, function(x) {
        match <- vals[grepl(paste0("^", x, "$"), nms)]
        res <- if(length(match) == 1) match else match[i]
        res
      }, USE.NAMES = FALSE)
      ##s6params(setNames(as.numeric(c(vls, p$n, p$epsilon_a, p$epsilon_r, p$A, p$eta_m)), c(parnms, "n", "epsilon_a", "epsilon_r", "A", "eta_m")))
      s6params(c(vls, p[c("n", "epsilon_a", "epsilon_r", "A", "eta_m")]))
    })
    Fmsy <- sapply(estpars, calcFmsy)
    SSBrel <- sapply(estpars, function(x) getParams(x)$B) / 
              mapply(function(p, fmsy) getParams(s6params(c(Fm = fmsy), base = p))$B, estpars, Fmsy, SIMPLIFY = TRUE)
    if(nyrs == 1) {
      estpars <- estpars[[1]]
      Fmsy <- Fmsy[[1]]
    }
    opt$convergence
  }, silent = !verbose)
  if(class(tryer) == "try-error") {
    return(tryer)
  }
  nw <- function(x) grepl(paste0("^", x, "$"), nms)
  structure(data.frame(Fm=vals[nw("Fm")], Fm_sd = sds[nw("Fm")],
                       Winf=rep(vals[nw("Winf")], nyrs), Winf_sd=rep(sds[nw("Winf")], nyrs),
                       Fmsy = Fmsy, FFmsy = vals[nw("Fm")]/Fmsy,
                       ##Fmsyest = vals[nw("Fmsy")], Fmsyest_sd = sds[nw("Fmsy")],
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
                       ssbrel = SSBrel,
                       row.names=inp@years),
            obj=obj, opt=opt, sdr = sdr, estpars=estpars)
}

