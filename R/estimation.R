
#' Title
#'
#' @param df 
#' @param n 
#' @param epsilon_a 
#' @param epsilon_r 
#' @param A 
#' @param eta_m 
#' @param a 
#' @param Winf 
#' @param sigma 
#' @param u 
#' @param sdloga 
#' @param winf.ubound 
#' @param Wfs 
#' @param verbose 
#' @param map 
#' @param random 
#' @param isSurvey 
#' @param eta_S 
#' @param usePois 
#' @param totalYield 
#' @param ... 
#'
#' @return
#' @export
#'
#' @examples
estimate <- function(df, n=0.75, epsilon_a=0.8, epsilon_r=0.1, A=4.47, 
                         eta_m=0.25, a=0.22, Winf = NULL, sigma=NULL, u = 10,
                         sdloga = 0.7, winf.ubound = 2, Wfs = NULL,
                         verbose = FALSE, map = list(loga = factor(NA)), 
                         random = c(), isSurvey = FALSE, eta_S = NULL, usePois = TRUE,
                         totalYield = NULL, ...) {
  if (is.null(df)) return(NULL)
  if (! require(TMB)) stop("TMB is not installed! Please install and try again.")
  isTS <- is(df, "list")
  if (isTS) {
    if (is.null(totalYield)) totalYield <- rep(0.01234567, length(df))
    length(totalYield) == length(df) || stop("Please provide the yield for all years")
    yrs <- names(df)
    df <- df2matrix(df)
    nyrs <- ncol(df)
    DLL <- "s6modelts"
  } else {
    if (is.null(totalYield)) totalYield <- 0.01234567
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
        eta_S <- which(apply(df, 1, function(x) sum(x) != 0))[[1]] * binsize / Winf
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
    if(is.null(u))  {
      u <- 10
    } else {
      map$logu  <- factor(NA)
    }
    if(isTS) {
      freq <- df
      nwc <- attr(df, "nwc")
      logFm <- rep(log(0.5), nyrs)
    } else {
      freq <- df$Freq
      nwc <- dim(df)[1]
      logFm <- log(0.5)
  }
    data <- list(binsize=binsize, nwc=nwc, freq=freq, n=n, epsilon_a=epsilon_a,
                 epsilon_r=epsilon_r, A=A, eta_m=eta_m, meanloga = log(a), 
                 sdloga = sdloga, isSurvey = as.integer(isSurvey),
                 usePois = as.integer(usePois), totalYield = totalYield)
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
    sdr <- sdreport(obj)
    parnms <- c("Fm","Winf","Wfs", "a", "eta_S")
    vals <- sdr$value
    nms <- names(vals)
    sds <- setNames(sdr$sd, nms)
    estpars <- sapply(seq(nyrs), function(i) {
      vls <- sapply(parnms, function(x) {
        match <- vals[grepl(paste0("^", x, "$"), nms)]
        if(length(match) == 1) match else match[i]
      })
      s6params(setNames(as.numeric(c(vls, n, epsilon_a, epsilon_r, A, eta_m)), c(parnms, "n", "epsilon_a", "epsilon_r", "A", "eta_m")))
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
                       row.names=yrs),
            obj=obj, opt=opt, sdr = sdr, estpars=estpars)
}

