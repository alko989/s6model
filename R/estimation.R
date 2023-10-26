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
#' sim <- simulateData3(params = parameters("Fm", 0.3, transformed=FALSE))
#'
#' ## Plotting the negative log likelihood for different values of Fm
#' Fm <- seq(0.1, 1, 0.05) 
#' Fm.transformed <- log(Fm / 0.25)
#' nll <- sapply(Fm.transformed, function(x) minimizeme(theta = x, data = sim$sample, names = "Fm"))
#' plot(Fm, nll, type="l") 
#'
#' ## Using optimise to estimate one parameter
#' est <- optimise(f = minimizeme, interval=c(0.01, 1), data=sim$sample, names="Fm")
#' est.Fm <- exp(est$minimum) * 0.25
#' abline(v=est.Fm)
#' mtext(paste("Estimated Fm = ", round(est.Fm, 2)), at=est.Fm)
#' }
#' 
#' @export minimizeme
#' @rdname minimizeme
minimizeme <- function(theta, data, names, fixed.names=c(), fixed.vals=c(), isSurvey=FALSE) {
  params <- parameters(c(names, fixed.names), c(theta, fixed.vals))  
  if(is(data, "data.frame")) {
    return(with(getParams(params,isSurvey),
                sum( - data$Freq * log(pdfN.approx(data$Weight)) )))
  } else if(is(data, "numeric")) {
    return(with(getParams(params,isSurvey), sum(-log(pdfN.approx(data)))))
  } else {
    stop("data appears not to be numeric vector or data.frame")
  }
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
estimate_TMB <- function(df, n=0.75, epsilon_a=0.8, epsilon_r=0.1, A=4.47, 
                         eta_m=0.25, a=0.22, Winf = NULL, sigma=NULL, u = NULL,
                         sdloga = 0.7, winf.ubound = 1.5, Wfs = NULL,
                         verbose=FALSE, map=list(), 
                         random=c(), isSurvey = FALSE, eta_S = NULL, usePois = TRUE,
                         totalYield = NULL, perturbStartingVals = FALSE, Fm = NULL,
                         Fmsy_range = FALSE, seed = NULL, it = 100, ...) {
  if (is.null(df)) return(NULL)
  if (! require(TMB)) stop("TMB is not installed! Please install and try again.") 
  isTS <- is(df, "list")
  if (isTS) {
    if(is.null(totalYield)) {totalYield <- rep(0.01234567, length(df))}
    length(totalYield) == length(df) || stop("Please provide the yield for all years")
    yrs <- names(df)
    df <- df2matrix(df)
    nyrs <- ncol(df)
    DLL <- "s6modelts"
  } else {
    if (is.null(totalYield)) {totalYield <- 0.01234567}
    yrs <- 1
    nyrs <- 1
    DLL <- "s6model"
  }
  tryer <- try({
    binsize <- attr(df,"binsize")
    if(is.null(Winf))  {
      if(isTS) {
        Winf <- (nrow(df) + 2) * binsize
      } else {
        Winf <- max(df$Weight) + 2 * binsize
      }
    } else {
      map$logWinf  <- factor(NA)
    }
    if(is.null(eta_S))  {
      if(isTS) {
        eta_S <- which(apply(df, 1, function(x) sum(x) != 0))[[1]] * binsize / Winf # Look at frequencies that are not 0, take the first, multiply by binsize to get actual size /Winf (look at definition) 
      
        } else {
        eta_S <- which(df$Freq != 0)[1] * binsize / Winf
      }
    } else {
      map$logeta_S  <- factor(NA)
    }
    if(! isSurvey) {
      map$logeta_S <- factor(NA)
    }    
    if(is.null(sigma) || is.na(sigma))  {
      if(isTS) {
        sigma <- if(usePois) colSums(df) else rep(0.0001, nyrs)
      } else {
        sigma <- if(usePois) sum(df$Freq) else 0.0001
      }
    } else {
      map$logSigma  <- rep(factor(NA), nyrs)
    }
    if(is.null(Wfs) || is.na(Wfs))  {
      if(isTS) {
        Wfs <- apply(df, 2, function(x) round((which.max(x) + which(x > 0)[1]) / 2) * binsize - binsize / 2)
      } else {
        Wfs <- df$Weight[round((which.max(df$Freq) + which(df$Freq > 0)[1]) / 2)] - binsize / 2## min(df$Weight[df$Freq > 0])
      }
     } else {
      map$logWfs  <- rep(factor(NA), nyrs)
    }
    # Calculate width of selectivity curve u:
    if(isTS) {
      if(!is.null(u))  {
        map$logu  <- factor(NA) # It's not estimated
        logu <- rep(log(u),nyrs) # Arbitrary initial value
      } else {
        logu  <- rep(log(10),nyrs)
      }
    } else {
      if(!is.null(u))  {
        map$logu  <- factor(NA) # It's not estimated
        logu <- log(u) # Arbitrary initial value
      } else {
        logu  <- log(10)
      }
    }
    # Calculate loga
    if(isTS) {
      if(!is.null(a))  {
        map$loga  <- factor(NA) 
        loga <- rep(log(a),nyrs)
      } else {
        loga  <- rep(log(0.22),nyrs)
      }
    } else {
      if(!is.null(a))  {
        map$loga  <- factor(NA) # It's not estimated
        loga <- log(a) # Arbitrary initial value
      } else {
        loga  <- log(0.22)
      }
    }
    # Calculate Fishing mortality Fm
    if(isTS) {
      freq <- df
      nwc <- attr(df, "nwc")
      if(!is.null(Fm)){
        map$logFm  <- factor(NA)
        logFm <- rep(log(Fm),nyrs)
      } else {
        logFm <- rep(log(0.5), nyrs)
      }
    } else {
      freq <- df$Freq
      nwc <- dim(df)[1]
      if(!is.null(Fm)){
        map$logFm  <- factor(NA)
        logFm <- log(Fm)
      } else {
        logFm <- log(0.5)
      }
    }
    data <- list(binsize=binsize, nwc=nwc, freq=freq, n=n, epsilon_a=epsilon_a,
                 epsilon_r=epsilon_r, A=A, eta_m=eta_m, meanloga = loga, 
                 sdloga = sdloga, isSurvey = as.integer(isSurvey),
                 usePois = as.integer(usePois), totalYield = totalYield)
    pars <- list(loga = loga, logFm = logFm, logWinf = log(Winf),
                 logWfs = log(Wfs), logSigma = log(sigma),logeta_S = log(eta_S), 
                 logu = logu)
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
    parnms <- c("Fm","Winf","Wfs", "a", "eta_S",'u','M') # Eta_f is included in Wfs
    vals <- sdr$value # This is not in log
    nms <- names(vals)
    sds <- setNames(sdr$sd, nms) # This is not in log
    # diag(sdr$cov.fixed) gives sds
    # sdr$cov.fixed
    estpars <- sapply(seq(nyrs), function(i) {
      vls <- sapply(parnms, function(x) {
        match <- vals[grepl(paste0("^", x, "$"), nms)]
        if(length(match) == 1) match else match[i]
      })
      parameters(c(parnms, "n", "epsilon_a", "epsilon_r", "A", "eta_m"),
                 as.numeric(c(vls, n, epsilon_a, epsilon_r, A, eta_m)),
                 transformed=FALSE)
    })
    
    if(Fmsy_range){
      
      # Set seed
      if(!is.null(seed))  {
        set.seed(seed)
      }
      
      # Draw from normally distributed parameters that are input to calcFmsy
      # Input params for calcFmsy:
      # Winf, n, epsilon_a, epsilon_r, A, eta_m, a, Wfs, u, F
      # Fm, Winf, a, Wfs and u are estimated -- somehow eta_S doesn't seem to be estimated?
      
      Fmsy.dist <- numeric(it)
      SSBrel.dist <- numeric(it)
      Bexplrel.dist <- numeric(it)
      
      # Clean up and remove for loop
      estpars.rand <- as.list(estpars[[1]])
      
      Fm <- rnorm(it, mean = as.numeric(estpars.rand[names(estpars.rand) %in% 'Fm']), sd = sds['Fm'])
      Winf <- rnorm(it, mean = as.numeric(estpars.rand[names(estpars.rand) %in% 'Winf']), sd = sds['Winf'])
      Wfs <- rnorm(it, mean = as.numeric(estpars.rand[names(estpars.rand) %in% 'Wfs']), sd = sds['Wfs'])
      a <- rnorm(it, mean = as.numeric(estpars.rand[names(estpars.rand) %in% 'a']), sd = sds['a'])
      u <- rnorm(it, mean = as.numeric(estpars.rand[names(estpars.rand) %in% 'u']), sd = sds['u'])
      
      for(i in 1:it){
        print(i)
        
        estpars.rand <- as.list(estpars[[1]])
        estpars.rand$Fm <- Fm[i]
        estpars.rand$Winf <- Winf[i]
        estpars.rand$Wfs <- Wfs[i]
        estpars.rand$a <- a[i]
        estpars.rand$u <- u[i]+
        estpars.rand <- parameters(c(names(estpars.rand)),
                                   as.numeric(c(estpars.rand)),
                                   transformed=FALSE)  
        
        # Fmsy[i] <- sapply(estpars.rand.1, calcFmsy)
        Fmsy.dist[i] <- calcFmsy(estpars.rand)
        # Bmsy: Biomass/Biomass at Fmsy
        SSBrel.dist[i] <- getParams(estpars.rand)$B / 
          getParams(parameters("Fm", Fmsy.dist[i], FALSE, base = estpars.rand))$B
        Bexplrel.dist[i] <-  getParams(estpars.rand)$Bexpl / 
          getParams(parameters("Fm", Fmsy.dist[i], FALSE, base = estpars.rand))$Bexpl
        
      }
      
      
      Fmsy <- sapply(estpars, calcFmsy)
      SSBrel <- sapply(estpars, function(x) getParams(x)$B) / 
        mapply(function(p, fmsy) getParams(parameters("Fm", fmsy, FALSE, base = p))$B, estpars, Fmsy, SIMPLIFY = TRUE)
      Bexplrel <- sapply(estpars, function(x) getParams(x)$Bexpl) / 
        mapply(function(p, fmsy) getParams(parameters("Fm", fmsy, FALSE, base = p))$Bexpl, estpars, Fmsy, SIMPLIFY = TRUE)
      if(nyrs == 1) {
        estpars <- estpars[[1]]
        Fmsy <- Fmsy[[1]]
      }
      
    } else {
      
      Fmsy <- sapply(estpars, calcFmsy)
      SSBrel <- sapply(estpars, function(x) getParams(x)$B) / 
        mapply(function(p, fmsy) getParams(parameters("Fm", fmsy, FALSE, base = p))$B, estpars, Fmsy, SIMPLIFY = TRUE)
      Bexplrel <- sapply(estpars, function(x) getParams(x)$Bexpl) / 
        mapply(function(p, fmsy) getParams(parameters("Fm", fmsy, FALSE, base = p))$Bexpl, estpars, Fmsy, SIMPLIFY = TRUE)
      if(nyrs == 1) {
        estpars <- estpars[[1]]
        Fmsy <- Fmsy[[1]]
      }
      
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
                       Wfs = vals[nw("Wfs")], Wfs_sd = sds[nw("Wfs")],
                       a = rep(vals[nw("a")], nyrs), a_sd = sds[nw("a")],
                       M = rep(vals[nw("M")], nyrs), M_sd = sds[nw("M")],
                       eta_S = vals[nw("eta_S")], eta_S_sd = sds[nw("eta_S")],
                       sigma = vals[nw("sigma")], sigma_sd = sds[nw("sigma")],
                       u = vals[nw("u")], u_sd = sds[nw("u")],
                       R = vals[nw("R")], R_sd = sds[nw("R")],
                       Rrel = vals[nw("Rrel")], Rrel_sd = sds[nw("Rrel")],
                       Rp = vals[nw("Rp")], Rp_sd = sds[nw("Rp")],
                       rmax = vals[nw("rmax")], rmax_sd = sds[nw("rmax")],
                       Y = vals[nw("Y")], Y_sd = sds[nw("Y")],
                       ssb = vals[nw("ssb")], ssb_sd = sds[nw("ssb")],
                       Bexpl = vals[nw("Bexpl")], Bexpl_sd = sds[nw("Bexpl")],
                       ssbrel = SSBrel, Bexplrel = Bexplrel, epsilon_a = epsilon_a, 
                       epsilon_r = epsilon_r, n = n, A = A, eta_m = eta_m,
                       row.names=yrs),
            obj=obj, opt=opt, sdr = sdr, estpars=estpars, log.sd = diag(sdr$cov.fixed))
}

# , Fmsy.dist = Fmsy.dist, SSBrel.dist = SSBrel.dist, Bexplrel.dist = Bexplrel.dist