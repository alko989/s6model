#' Returns a list with all parameters and functions
#' 
#' Given a s6params object, it returns as a list all relevant parameters and
#' corresponding functions like pdf and cdf.
#' 
#' 
#' @param p A s6params object
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
#' @seealso \code{\link{s6params}}
#' @keywords misc
#' @examples
#' 
#'  getParams()
#' 
#' @export
getParams <- function(p = s6params(),  FF=NULL, calcBRPs=FALSE, isSurvey=FALSE, 
                      optim.fmsy=FALSE, optim.fmsyr = FALSE, optim.Rrel =FALSE) {
  if( ! is(p, "s6params"))
    stop("Wrong input argument in getParams. Use the 's6params' class instead.")  
  with(as.list(p), {
    if (optim.fmsy | optim.fmsyr | optim.Rrel)
      Fm <- FF
    eta_F <- Wfs/Winf
    if (! isTRUE(all.equal(p@logeta_F, log(eta_F)))) warning("Wfs and eta_F do not match! Wfs was used and eta_F was returned correctly.")
    p@logeta_F <- log(eta_F)
    
    w_r <- w_egg <- 0.001
    Delta <- (log(Winf) - log(w_r)) / (M - 1)
    w <- exp(log(w_r) + (1:M - 1) * Delta)
    
    delta <- diff(w)
    
    psi_F <- (1 + (w / Wfs) ^ -u ) ^ -1
    psi_S <- (1 + (w / (eta_S * Winf)) ^ -u) ^ -1
    psi_m <- (1 + (w / (eta_m * Winf))^ -10 )^-1
    m <- a * A * w^(n - 1) +  Fm * psi_F
    m[M] <- Inf
    g <- A * w ^ n * (1 - (w / Winf) ^ (1 - n) * (epsilon_a + (1 - epsilon_a) * psi_m))
    
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
      pdfN.approx <- approxfun(w, pdfN, yleft = 0, yright = 0)## .Machine$double.xmin)
      cdf <- approxfun(w, cumsum(pdfN * c(0, delta)), yleft = 0, yright = 1)
    }
    B <- sum((psi_m  * N * w)[-M] * delta)
    Rrel <- 1 - (Winf^(1-n) * w_egg/(epsilon_r * (1 - epsilon_a) * A * B))## * (w_r/w_egg)^(a-1)
    Rp <- epsilon_r * (1 - epsilon_a) * A * Winf ^ (n-1) / w_egg * B
    Y <- Fm * Rrel * sum((psi_F * N * w)[-M] * delta)
    wF <- which.min(abs(Wfs - w)) ## Find the closest weight to Wfs
    
    YR <- Y * 1 / (N[wF] * Rrel * g[wF])# Fm * sum((psi_F * N * w)[-M] * delta)
    
    if(calcBRPs) {
      Fmsy <- optimise(f = getParams, interval=c(0,10), maximum=TRUE, 
                       p = p, optim.fmsy = TRUE)$maximum
      FoverFmsy <- Fm / Fmsy
      ##Fmsyr <- optimise(f=getParams, interval=c(0,2), maximum=TRUE, p=p, optim.fmsyr=TRUE)$maximum
      ##FoverFmsyr <- Fm/Fmsyr
      Fcrash <- try(uniroot(f = getParams, interval=c(1e-25,10),
                            p = p, optim.fmsy = TRUE)$root, silent = TRUE)
      Flim <-try(uniroot(f = getParams, interval = c(1e-25, 10),
                         p = p,  optim.Rrel = TRUE)$root, silent = TRUE)
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
  })
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
##' @return \code{s6params}
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
      res <- s6params(setNames(par.vals, parameter.names))
      if( getParams(res)$Rrel >=  Rrel.gt & getParams(res, calcBRPs=TRUE)$Fmsy >= Fmsy.gt) return(res)
    }
  }

##' Random parameters with fixed Winf
##'
##' Shorthand function for a specific asymptotic weight
##' @param winf Numeric. Asymptotic weight
##' @return \code{s6params}
##' @author alko
##' @rdname getRandomParameters
##' @export
getRandomParameters.fixedWinf <- function(winf, Rrel.gt=-Inf, Fmsy.gt=0) {
  parameter.names <- c("A", "n" ,"eta_m","eta_F", "a" ,"Fm","Winf","epsilon_a", "epsilon_r")
  parameter.value <- c(4.5,0.75 , 0.25  ,  0.05 , 0.35,0.25,  winf ,    0.8    ,     0.1    )
  getRandomParameters(parameter.names, parameter.value, Rrel.gt=Rrel.gt, Fmsy.gt=Fmsy.gt)
}


##' Calculates the Fmsy reference point
##'
##' @param params An object of class \code{s6params} 
##' @return numeric F that leads to MSY
##' @author alko
##' @export
calcFmsy <- function(params=NULL) {
  if(!require(TMB)) stop("TMB package not installed.")
  if(is.null(params)) return (NULL)
  if(is(params, "s6params")) {
    params <- as.list(params)
  }
  if( ! is(params,"list"))
    stop("params is of class ", class(params))
  def <- list(n=0.75, epsilon_a=0.8, epsilon_r=0.1, A=4.47, eta_m=0.25, a=0.27, M=1000, u = 10)
  def <- replace(def, names(params), unlist(params))
  obj <- MakeADFun(def, list(logF = log(0.2)), DLL="calcFmsy")
  obj$env$tracemgc <- FALSE
  obj$env$silent <- TRUE
  newtonOption(obj=obj,trace=0); config(trace.optimize = 0, DLL="calcFmsy")
  opt <- try(do.call("optim", obj))
  res <- try(sdreport(obj)$val)
  if(is(res, "try-error"))
    return(NULL)
  names(res) <- NULL
  res
}

##' Change the bin size of a weight frequency data.frame
##' 
##' 
##' @param df data.frame with colums 'Weight' and 'Freq'
##' @param binsize numeric, the bin size in grams
##' @param weight.col The weight column name
##' @param freq.col The frequency column name
##' @param keepZeros logical, if FALSE the bins with zero observations are removed
##'
##' @export
##' 
##' @examples 
##' ## Simulate a data set with bin size equal to 100 gr
##' dat <- simulateData3(binsize=100)
##' 
##' ## Change the bin size to 200 gr
##' dat <- changeBinsize(dat$df, binsize = 200)
##' @return data.frame with weight frequencies binned with bin size \code{binsize}
changeBinsize <- function(df, binsize = 10, keepZeros = TRUE, weight.col = "Weight", freq.col = "Freq") {
  if(is(df, "list")) { 
    return(lapply (df, changeBinsize2, binsize = binsize, keepZeros = keepZeros, 
                   weight.col = weight.col, freq.col = freq.col)) 
  }
  if(dim(df)[1] == 0 ) return( structure(data.frame(Weight = numeric(0), Freq = integer(0)), binsize = binsize))
  if(require(dplyr)) {
    mx <- max(df[weight.col], na.rm = TRUE)
    res <- df %>% 
      mutate_(Weight = weight.col) %>% 
      mutate(Weight = Weight - (Weight - 1) %% binsize + binsize / 2 - 1) %>%
      group_by(Weight) %>%  
      summarise_(Freq = ~as.integer(round(sum(Freq)))) %>%
      full_join(data.frame(Weight = seq(binsize / 2, mx + binsize, binsize)), by = "Weight") %>%
      mutate(Freq = as.integer(ifelse(is.na(Freq), 0, Freq)))
  } else {
    cuts <- seq(0, max(df[weight.col], na.rm = TRUE) + binsize, binsize)
    labs <- head(cuts + binsize/2, -1)
    res <- data.frame(Weight = c(cut(df[[weight.col]], breaks = cuts, labels = labs), labs), 
                      Freq = c(df[[freq.col]], rep(0, length(labs))),
                      stringsAsFactors = TRUE)
    res <- aggregate(Freq ~ Weight, data = res, FUN = sum)
    res$Weight <- as.numeric(res$Weight)
    res$Freq <- as.integer(round(as.numeric(res$Freq), digits = 0))
    rl <- rle(res$Freq)
    res <- head(res, nrow(res) - if(tail(rl$values, 1) == 0) tail(rl$lengths, 1) else 0)  
  }
  if (! keepZeros) {
    res <- res[df$Freq > 0, ]
  }
  res <- structure(res, binsize = binsize)
  res
}

##' @export
##' @title Find newest file in a folder matching a pattern
##' @description Given a folder and a pattern rerurns the file that was modified latest.
##'
##' @param path character, relative or absolute path
##' @param pattern the pattern to search for, used by \code{link{dir}}
##'
##' @seealso \code{\link{dir}}
findLatest <- function(path = ".", pattern = "") {
  files <- dir(path, pattern, ignore.case = TRUE, full.names = TRUE)
  if(length(files) == 0) return(NULL)
  if(length(files) == 1) return(files)
  df <- lapply(files, function(f){
    data.frame(fn = f, mtime = file.info(f)$mtime, stringsAsFactors = FALSE)
  })
  df <- do.call(rbind.data.frame, df)
  latest <- df[order(df$mtime, decreasing = TRUE), ][1,1]  
  latest
}
