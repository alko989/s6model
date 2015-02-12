#' Returns a list with all parameters and functions
#' 
#' Given a Parameters object, it returns as a list all relevant parameters and
#' corresponding functions like pdf and cdf.
#' 
#' 
#' @param p A Parameters object
#' @param FF Numeric. Fishing mortality. This argument is ignored if all optim.* are FALSE 
#' @param calcBRPs Boolean. If true, calculates biological reference points.
#' @param isSurvey boolean, if TRUE a survey selectivity is used for the pdf
#' @param optim.fmsy Logical.
#' @param optim.fmsyr Logical.
#' @param optim.Rrel Logical
#' @return \itemize{
#' \item If optim.* are all FALSE, an invisible list, containing all model parameters, the cdf, the pdf
#' functions. spawning stock biomass (SSB), yield. If calcBRPs is TRUE also Fmsy.
#' \item If optim.fmsy is TRUE, only the yield per Rmax is returned
#' \item If optim.fmsyr is TRUE, only the yield per recruit is returned
#' \item If optim.Rrel is TRUE, only the R/Rmax - 0.5 is returned
#' }
#' @author alko
#' @seealso \code{\link{parameters}}
#' @keywords misc
#' @examples
#' 
#'  getParams()
#' 
#' @export
getParams <- function(p = new("Parameters"),  FF=NULL, calcBRPs=FALSE, isSurvey=FALSE, 
                      optim.fmsy=FALSE, optim.fmsyr = FALSE, optim.Rrel =FALSE) {
  if( ! is(p,"Parameters"))
    stop("Wrong input argument in getParams. Use the Parameters class instead.")  
  Winf <- exp(p@logWinf) * p@scaleWinf
  Fm <- exp(p@logFm) * p@scaleFm
  if(optim.fmsy | optim.fmsyr | optim.Rrel)
    Fm <- FF
  A <- exp(p@logA) * p@scaleA
  n <- exp(p@logn) *p@scalen
  eta_m <- exp(p@logeta_m) * p@scaleeta_m
  eta_S <- exp(p@logeta_S) * p@scaleeta_S
  a <- exp(p@loga) *  p@scalea
  epsilon_a <- exp(p@logepsilon_a) * p@scaleepsilon_a
  epsilon_r <- exp(p@logepsilon_r) * p@scaleepsilon_r
  Wfs <- exp(p@logWfs) * p@scaleWfs
  eta_F <- Wfs/Winf
  if(!(isTRUE(all.equal(p@logeta_F, log(eta_F / p@scaleeta_F))))) warning("Wfs and eta_F do not match! Wfs was used and eta_F was returned correctly")
  p@logeta_F <- log(eta_F / p@scaleeta_F)
  u <-exp(p@logu) * p@scaleu
  M <- p@M
  
  w_r <- w_egg<- 0.001
  Delta <- (log(Winf) - log(w_r)) / (M - 1)
  w <- exp(log(w_r) + (1:M - 1) * Delta)
  
  delta <- diff(w)
  
  psi_F <- (1 + (w / Wfs)^-u )^-1
  
  psi_S <- (1 + (w / (eta_S * Winf))^-u )^-1
  
  psi_m <- (1 + (w / (eta_m * Winf))^-10 )^-1
  m <- a * A * w^(n - 1) +  Fm * psi_F
  m[M] <- Inf
  g <- A*w^n *(1 - (w/Winf)^(1-n) * (epsilon_a + (1 - epsilon_a) * psi_m))
  
  N <- exp(- cumsum((m / g)[-1] * delta)) / g[-1]
  N <- c(1/g[1], N)
  N[M] <- 0
  if(isSurvey)
  {
    fishing <- psi_S * N
  }
  else
  {
    fishing <- psi_F * N
  }
  if( ! (optim.fmsy | optim.fmsyr | optim.Rrel)) {    
    pdfN <-  fishing / sum(fishing * c(delta, 0))
    pdfN.approx <- approxfun(w, pdfN, yleft=0, yright=.Machine$double.xmin)
    cdf <-approxfun(w, cumsum(pdfN * c(delta,0)), yleft=0, yright=1)
  }
  B <- sum((psi_m  * N * w)[-M] * delta)
  Rrel <- 1 - (Winf^(1-n) * w_egg/(epsilon_r * (1 - epsilon_a) * A * B))## * (w_r/w_egg)^(a-1)
  
  Y <- Fm * Rrel * sum((psi_F * N * w)[-M] * delta)
  YR <- Fm * sum((psi_F * N * w)[-M] * delta)
  
  if(calcBRPs) {
    Fmsy <- optimise(f=getParams, interval=c(0,10), maximum=TRUE, p=p, optim.fmsy=TRUE)$maximum
    FoverFmsy <- Fm/Fmsy
    ##Fmsyr <- optimise(f=getParams, interval=c(0,2), maximum=TRUE, p=p, optim.fmsyr=TRUE)$maximum
    ##FoverFmsyr <- Fm/Fmsyr
    Fcrash <- try(uniroot(f=getParams, interval=c(1e-25,10), p=p, optim.fmsy=TRUE)$root, silent=TRUE)
    Flim <-try(uniroot(f=getParams, interval=c(1e-25,10), p=p,  optim.Rrel = TRUE)$root, silent=TRUE)
  }   
  vb.M <- a * A * Winf^(n-1)*eta_m^(n-1)
  vb.K <- A * Winf^(n-1) / 3
  vb.MK <- vb.M / vb.K
  if(optim.fmsy)
    return(Y)
  if(optim.fmsyr)
    return(YR)
  if(optim.Rrel)
    return(Rrel - 0.5)
  res <- as.list(environment())
  attr(res, "version") <- getVersion()
  return(invisible(res))
}

##' Simulates catch-at-weight data using the s6model
##'
##' Simulate a data set of individual catch weight or 
##' @param samplesize Integer. Number of simulated individuals
##' @param params Object of class \code{Parameters}. Model parameters
##' @param binsize Numeric. Weight class width. If \code{retDF} is FALSE, \code{binsize} is ignored
##' @param keepZeros Logical. If TRUE the resulting data.frame includes weight classes with zero individuals. Otherwise these weight classes are dropped. Ignored if \code{retDF} is FALSE
##' @param retDF Logical. If TRUE a data.frame is returned with columns Weight (weight classes) and Freq (numbers per weight class)
##' @param ... Extra named arguments passed to \code{getParams}
##' @return Invisible list with
##' \itemize{
##'   \item `sample` containing weights of individuals,
##'   \item `parameters` the parameters as a \code{Parameters} object
##'   \item `Fmsy` the Fmsy of the given parameters
##'   \item `df` A data.frame with columns Weight and Freq with number per weight class. This is returned only if \code{retDF} is TRUE
##' }
##' @author alko
##' @export
simulateData3 <- function(samplesize= 1000, params = parameters(), binsize = 5, keepZeros=TRUE, retDF=TRUE, ...)
{
  applyfun <- if(require(parallel)) mclapply else sapply
  sam <- c()
  with(getParams(params, ...), {
    sam <<- simplify2array(applyfun(runif(samplesize), function(u) uniroot(function(x) {cdf(x) - u}, c(w_r, Winf))$root))
  })
  res <- list(sample = sam, parameters = params, Fmsy = getParams(params,calcBRPs=TRUE)$Fmsy)
  if(retDF) res$df <- sample2df(sam, binsize, keepZeros=keepZeros)
  attr(res, "version") <- getVersion()
  return(invisible(res))  
}

##' Convert a vector sample to data.frame with counts per weight class
##'
##' Takes a vector containing individual catch weights and returns a data.frame with numbers per weight class
##' @param sam Numeric vector. Individual weights
##' @param binsize Numeric. Weight class width
##' @param keepZeros Logical. If TRUE the resulting data.frame includes weight classes with zero individuals. Otherwise these weight classes are dropped.
##' @return A data.frame with columns Weight and Freq with number per weight class.
##' @author alko
##' @export
sample2df <- function(sam, binsize, keepZeros=TRUE) {
  df <- as.data.frame(table(cut(sam, seq(0,max(sam) + binsize, binsize),
                                labels=seq(binsize/2, max(sam) + binsize/2, binsize) )), stringsAsFactors=FALSE)
  names(df) <- c("Weight","Freq")
  if (! keepZeros) df <- df[df$Freq > 0, ]
  df$Weight <- as.numeric(df$Weight)
  attr(df, "binsize") <- binsize
  df
}

rparam <- function(value, range.sd, range.cv, lb= -Inf, ub = Inf, unif=FALSE)
{
  if(unif)
    return(runif(1,lb, ub))
  if(range.sd == 0) {
    if(range.cv == 0) {
      return(value)
    } else {
      res <- rnorm(1,value, value*range.cv)
      while( (res < lb) || (res > ub) ) {
        res <- rnorm(1,value, value*range.cv)
      }
      return(res)
    }
  } else {
    res <- exp(rnorm(1, log(value), range.sd))
    while((res < lb) || (res > ub) ) {
      res <- rlnorm(1, log(value), range.sd)
    }
    return(res)    
  }
}

##' Random model parameters
##'
##' Random model parameters with constraints in Fmsy and relative recruitment (Rrel = R/Rmax)
##' @param parameter.names Character vector. Names of the parameters
##' @param parameter.value Numeric vector. Parameter mean values
##' @param parameter.sd Numeric vector. Standard deviation of the parameters
##' @param parameter.cv Numeric vector. Coefficient of variation of log transformed parameters
##' @param parameter.lbound Numeric vector. Lower bound of the distributions
##' @param parameter.ubound Numeric vector. Upper bound of the distributions
##' @param parameter.unif Logical vector. Use uniformly distributed parameters
##' @param Rrel.gt Numeric. Relative recruitment constraint. It allows parameters that lead
##' to Rrel at least equal to Rrel.gt
##' @param Fmsy.gt Numeric. Fmsy constraint. It allows parameters that lead
##' to Fmsy at least equal to Fmsy.gt
##' @return \code{Parameters}
##' @author alko
##' @export
getRandomParameters <-
  function(parameter.names=c("A", "n" ,"eta_m","eta_F", "a" ,"Fm","Winf","epsilon_a", "epsilon_r"),
           parameter.value =c(4.5,0.75 , 0.25  ,  0.05 , 0.35,0.25,  1e4 ,    0.8    ,     0.1    ),
           parameter.sd    =c( 0 ,  0  ,   0   ,  0.5  ,  0  , 0  ,   0  ,     0     ,      0     ),
           parameter.cv    =c(0.5,  0  ,  0.3  ,   0   , 0.5 , 0  ,   0  ,    0.1    ,     0.5    ),
           parameter.lbound=c( 0 ,  0  ,  0.01 ,   0   ,0.01 , 0  ,   0  ,     0     ,      0     ),
           parameter.ubound=c(Inf,  0  ,  Inf  ,  Inf  , Inf , 0  ,   0  ,    0.99   ,      Inf   ),
           parameter.unif  =c(F,F,F,F,F,F,F,F,F), Rrel.gt=-Inf, Fmsy.gt=0) {
    while(TRUE) {
      par.vals <- sapply(seq(along.with=parameter.value),
                         function(x) rparam(parameter.value[x], 
                                            parameter.sd[x], 
                                            parameter.cv[x],
                                            parameter.lbound[x],
                                            parameter.ubound[x],
                                            parameter.unif[x]))
      res <- parameters(parameter.names, par.vals, FALSE)
      if( getParams(res)$Rrel >=  Rrel.gt & getParams(res, calcBRPs=T)$Fmsy >= Fmsy.gt) return(res)
    }
  }

##' Random parameters with fixed Winf
##'
##' Shorthand function for a specific asymptotic weight
##' @param winf Numeric. Asymptotic weight
##' @return \code{Parameters}
##' @author alko
##' @rdname getRandomParameters
##' @export
getRandomParameters.fixedWinf <- function(winf, Rrel.gt=-Inf, Fmsy.gt=0) {
  parameter.names <- c("A", "n" ,"eta_m","eta_F", "a" ,"Fm","Winf","epsilon_a", "epsilon_r")
  parameter.value <- c(4.5,0.75 , 0.25  ,  0.05 , 0.35,0.25,  winf ,    0.8    ,     0.1    )
  getRandomParameters(parameter.names, parameter.value,, Rrel.gt=Rrel.gt, Fmsy.gt=Fmsy.gt)
}

##' Multicore lapply function with progress bar
##'
##' Wrapper around parallel::mclapply function with text progress bar
##' @param X  a vector (atomic or list) or an expressions vector.  Other
##' objects (including classed objects) will be coerced by `as.list`.
##' @param FUN the function to be applied to (`mclapply`) each element of `X` or (`mcmapply`) in parallel to `...`.
##' @param ... Optional arguments to FUN
##' @param progressbar Logical. If TRUE a text progress bar is shown.
##' @return Result from mclapply
##' @note If mclapply is not available, lapply is used instead.
##' @author alko
##' @export
tmclapply <- function(X, FUN, ..., progressbar=TRUE){
  aplfun <- if(require(parallel)) mclapply else lapply
  start <- Sys.time()
  if(progressbar)
    pb <- txtProgressBar(min = 0, max = 100, style=3)
  results <- local({
    f <- fifo(tempfile(), open="w+b", blocking=TRUE)
    if (inherits(parallel:::mcfork(), "masterProcess")) {
      progress <- 0.0
      while(progress < 1 && !isIncomplete(f)) {
        msg <- readBin(f, "double")
        progress <- progress + as.numeric(msg)
        if(progressbar)
          setTxtProgressBar(pb, progress * 100)
        tt <- (1 - progress)*(difftime(Sys.time(), start, units="mins"))/ progress
        cat(" ETC:", as.integer(tt), "min(s) and", round((tt - as.integer(tt)) * 60, 0) ,"secs")
        if( ! progressbar) cat("\r")
      } 
      parallel:::mcexit()
    }
    res <- aplfun(X, function(x) {
      rr <- FUN(x)
      writeBin(1/length(X), f)
      rr
    })
    close(f)
    if(progressbar) {
      setTxtProgressBar(pb,100)
      close(pb)
    }
    res
  })
  cat(difftime(Sys.time(), start, units="mins"), "mins\n")
  results
} 

##' Calculates the Fmsy reference point
##'
##' @param params An object of class \code{Parameters} 
##' @return numeric F that leads to MSY
##' @author alko
##' @export
calcFmsy <- function(params=NULL) {
  if(!require(TMB)) stop("TMB package not installed.")
  if(is.null(params)) return (NULL)
  if(is(params, "Parameters")) {
    params <- as.list(params)
  }
  if( ! is(params,"list"))
    stop("params is of class ", class(params))
  def <- list(n=0.75, epsilon_a=0.8, epsilon_r=0.1, A=4.47, eta_m=0.25, a=0.35, M=1000)
  def <- replace(def, names(params), unlist(params))
  obj <- MakeADFun(def, list(logF = log(0.2)), DLL="calcFmsy")
  obj$env$tracemgc <- FALSE
  obj$env$silent <- TRUE; newtonOption(trace=0); config(trace.optimize = 0,DLL="calcFmsy")
  opt <- try(do.call("optim", obj))
  res <- try(sdreport(obj)$val)
  if(is(res, "try-error"))
    return(NULL)
  names(res) <- NULL
  res
}
## simplify2dataframe <- function(dd) {
##   data.frame(t(simplify2array(dd)))
## }

##' @title Change the bin size of a weight frequency data.frame
##' 
##' @param df data.frame with colums 'Weight' and 'Freq'
##' @param binsize numeric, the bin size in grams
##' @param keepZeros logical, if TRUE the bins with zero observation are not removed
##' @export
##' 
##' @examples 
##' ## Simulate a data set with bin size equal to 100 gr
##' dat <- simulateData3(binsize=100)
##' 
##' ## Change the bin size to 200 gr
##' dat <- changeBinsize(dat$df, binsize = 200)
##' @return data.frame with weight frequencies binned with bin size \code{binsize}
changeBinsize <- function(df, binsize = 10, keepZeros = TRUE) {
  sample2df(rep(df$Weight, df$Freq), binsize, keepZeros=keepZeros)
}

##' @export
changeBinsize2 <- function(df, binsize = 10, keepZeros = TRUE, weight.col = "Weight", freq.col = "Freq", verbose = options()$verbose) {
  cuts <- seq(0, max(df[weight.col], na.rm = TRUE) + binsize, binsize)
  labs <- head(cuts + binsize/2, -1)
  res <- data.frame(Weight = cut(df[[weight.col]], breaks = cuts, labels = labs), Freq = df[[freq.col]])
  res <- rbind(res, data.frame(Weight = labs, Freq = 0))
  res <- aggregate(Freq ~ Weight, data = res, FUN = sum )
  res$Weight <- as.numeric(as.character(res$Weight))
  res$Freq <- as.numeric(as.character(res$Freq))
  if (! keepZeros) {
    res <- res[df$Freq > 0, ]
  }
  res <- structure(res, binsize = binsize)
  if(verbose) cat("Done")
  res
}

##' @export
##' @title Find newest file in folder matching pattern
##' @description Given a folder and a pattern, rerurns the file that was modified latest.
##' @param path character, relative or absolute path
##' @param patterm an optional \code{\link{regex}}. Only file names which match the regular expression will be returned. Passed to \code{\link{dir}}
##' @seealso \code{dir}
findLatest <- function(path = ".", pattern = "") {
  files <- dir(path, pattern, ignore.case = TRUE, full.names = TRUE)
  
  df <- lapply(files, function(f){
    data.frame(fn = f, mtime = file.info(f)$mtime, stringsAsFactors = FALSE)
  })
  df <- do.call(rbind.data.frame, df)
  latest <- df[order(df$mtime, decreasing = TRUE), ][1,1]  
  latest
}