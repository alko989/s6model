## Function for calculation of the steady-state single-species size-spectrum
## Author: Alexandros Kokkalis
## 
## Work in constant progress
## Changelog
## * 2012-6-22 Made function getParams that calculates the size-spectrum,
##             growth, mortality, SSB and yield
## * 2012-6-23 Sampling function
## * 2012-6-26 Added pdf and cdf functions
##  *** The changelog is moved to an extra file ***

getParams <- function(p = new("Parameters"), FF=NULL,isSurvey=FALSE, calcBRPs=FALSE, optim.fmsy=FALSE, optim.fmsyr = FALSE, optim.Rrel =FALSE)
  {
    if(class(p) != "Parameters")
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
        fishing <- psi_F * Fm * N
      }
    if( ! (optim.fmsy | optim.fmsyr | optim.Rrel)) {    
      pdfN <-  fishing / sum(fishing * c(delta, 0))
      pdfN.approx <- approxfun(w, pdfN, yleft=0, yright=0)
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
    else if(optim.fmsyr)
        return(YR)
    else if(optim.Rrel)
        return(Rrel - 0.5)
    else
        return(invisible(as.list(environment())))
  }

getParams2 <- function(p = new("Parameters"), isSurveY=FALSE)
    {
    if(class(p) != "Parameters")
      stop("Wrong input argument in getParams. Use the Parameters class instead.")
    Winf <- exp(p@logWinf) * p@scaleWinf
    Fm <- exp(p@logFm) * p@scaleFm
    A <- exp(p@logA) * p@scaleA
    n <- exp(p@logn) *p@scalen
    eta_F <- exp(p@logeta_F) * p@scaleeta_F
    eta_m <- exp(p@logeta_m) * p@scaleeta_m
    a <- exp(p@loga) *  p@scalea
    epsilon_a <- exp(p@logepsilon_a) * p@scaleepsilon_a
    epsilon_r <- exp(p@logepsilon_r) * p@scaleepsilon_r
    M <- 300
   
    w_r <- w_egg<- 0.001
    Delta <- (log(Winf) - log(w_r)) / (M - 1)
    w <- exp(log(w_r) + (1:M - 1) * Delta)
    
    wrel <- w/Winf
    delta <- diff(c(w, Winf))
    
    u <- 10
    psi_F <- (1 + (w / (eta_F * Winf))^-u )^-1
    ## psi_F.func <- function(ww) (1 + (ww / (eta_F * Winf))^-u)^-1
    psi_m <- (1 + (w / (eta_m * Winf))^-u )^-1
    ## psi_m.func <- function(ww) (1 + (ww / (eta_m * Winf))^-u)^-1
    m <- a * A * w^(n - 1) +  Fm * psi_F
    m[M] <- Inf
    ## m.func <- function(ww) a*h*ww^(n-1) + Fm * psi_F.func(ww)
    g <- A*w^n *(1 - (w/Winf)^(1-n) * (epsilon_a + (1 - epsilon_a) * psi_m))
    ##g.func <- function(ww) h*ww^n*(1-(ww/Winf)^(1-n) * (epsilon_a + (1 - epsilon_a) * psi_m.func(ww)))
    N <- cumprod(exp(Delta)^(-m/g*w))/(g*w_r)
    ## N.approx <- approxfun(w, N, yleft = 0, yright = 0) 
    ## if(require(multicore)) {
    ##  N.func <- function(x) unlist(mclapply(x, function(ww) exp(-integrate(function(x) m.func(x)/g.func(x), lower = w_r, upper = ww)$value) / g.func(ww)))
    ## } else {
    ##  N.func <- function(x) unlist(lapply(x, function(ww) exp(-integrate(function(x) m.func(x)/g.func(x), lower = w_r, upper = ww)$value) / g.func(ww))) 
    ## }
    ## N1 <- N[1]
    ## N <- N / N1
    N[M] <- 0
    if(isSurvey)
      {
        psiSurvey <- 1
        fishing <- psiSurvey * N
        pdfN <- fishing / sum(fishing * delta)
      }
    else
      {
        fishing <- psi_F * Fm * N
        pdfN <-  fishing / sum(fishing * delta)
      }
    pdfN.approx <- approxfun(w, pdfN, yleft=0, yright=0)
    ## pdfN.func <- function(ww) fishing.func(ww) / fishsum
    ## cdf <-approxfun(w, cumsum(pdfN), yleft=0, yright=0)
    cdf <-approxfun(w, cumsum(pdfN * delta), yleft=0, yright=1)
    B <- sum(psi_m  * N * w * delta)
    ## Rrel <-(1 - (Winf^(1-n) * w_r^(1+n))/(epsilon_r*B))/ h * w_r^(n+1)
    Rrel <- 1 - Winf^(1-n)/(epsilon_r * (1 - epsilon_a) * A * B) * (w_r/w_egg)^(a-1)
    ##Y <- Fm * Rrel* sum(psi_F * N  * w * delta) 
    Y <- Fm * Rrel * sum(psi_F * N * w * delta)
    #N0 <- exp(- cumsum((a * A * w^(n - 1)) / g * delta)) / g
    ##N0 <- N0 / N0[1]
    #norm = A*(eta_F * Winf)^(n+1) * (1 - (eta_F)^(1 - n))* N0[which(w >= (eta_F * Winf))][1]
    #YpR = Fm * sum(psi_F * N  * w * delta)  / norm
    Fmsy <- optimise(f=getY2, interval=c(0,10), maximum=TRUE, p=p)$maximum
    FoverFmsy <- Fm/Fmsy
    return(invisible(as.list(environment())))
  }




simulateData <- function(samplesize= 1000, params = parameters())
{
  sam <- c()
  with(getParams(params), {
    cdf.sim <- function(F, ...) {
      res <- NULL
      i <- 0
      while(is.null(res))
        {
          U <- runif(1)
          res <- tryCatch(uniroot(function(x) {cdf(x) - U}, c(w_r, Winf))$root,
                          error=function(e)
                          {
                            i <- i + 1
                            if(i > 999)
                              warning(" ** No sample found after 1000 tries **")
                            warning(paste("Error message: ",  e$message))
                            return(NULL)
                          })
        }
      return(res)
    }
    
    for(i in 1:samplesize)
      sam <<-  c(sam, cdf.sim(cdf))

  })
  return(sam)
}

simulateData2 <- function(samplesize= 1000, params = parameters(), ...)
{
  sam <- c()
  with(getParams(params, ...), {
    cdf.sim <- function() {
      res <- NULL
      i <- 0
      while(is.null(res))
      {
        U <- runif(1)
        res <- uniroot(function(x) {cdf(x) - U}, c(w_r, Winf))$root
      }
      return(res)
    }
    sam <<- simplify2array(mclapply(1:samplesize, function(x) cdf.sim()))
  })
  table <- as.data.frame(table(cut(sam, seq(0,20000, 50), labels=seq(25, 20000 - 25, 50) )))
  names(table) <- c("Weight","Freq")
  table <- table[table$Freq > 0, ]
  return(list(sample = sam, parameters = params, Fmsy = getParams(params,calcBRPs=TRUE)$Fmsy, table= table))
}

simulateData3 <- function(samplesize= 1000, params = parameters(), wcw = 5, keepZeros=TRUE, retDF=TRUE, ...)
{
  applyfun <- if(require(multicore)) mclapply else sapply
  sam <- c()
  with(getParams(params, ...), {
    sam <<- simplify2array(applyfun(runif(samplesize), function(u) uniroot(function(x) {cdf(x) - u}, c(w_r, Winf))$root))
  })
  res <- list(sample = sam, parameters = params, Fmsy = getParams(params,calcBRPs=TRUE)$Fmsy)
  if(retDF) res$df <- sample2df(sam, wcw, keepZeros=keepZeros)
  return(res)  
}

sample2df <- function(sam, wcw, keepZeros=TRUE) {
  df <- as.data.frame(table(cut(sam, seq(0,max(sam) + wcw, wcw),
                                labels=seq(wcw/2, max(sam) + wcw/2, wcw) )), stringsAsFactors=FALSE)
  names(df) <- c("Weight","Freq")
  if (! keepZeros) df <- df[df$Freq > 0, ]
  df$Weight <- as.numeric(df$Weight)
  attr(df, "wcw") <- wcw
  df
}

getY <- function(F, p)
  {
    Winf <- exp(p@logWinf) * p@scaleWinf
    A <- exp(p@logA) * p@scaleA
    n <- exp(p@logn) *p@scalen
    Wfs <- exp(p@logWfs) * p@scaleWfs
    eta_m <- exp(p@logeta_m) * p@scaleeta_m
    a <- exp(p@loga) *  p@scalea
    epsilon_a <- exp(p@logepsilon_a) * p@scaleepsilon_a
    epsilon_r <- exp(p@logepsilon_r) * p@scaleepsilon_r
    M <- p@M
    w_r <- w_egg <- 0.001
    Delta <- (log(Winf) - log(w_r)) / (M - 1)
    w <- exp(log(w_r) + (1:M - 1) * Delta)
    rm(Delta)
    wrel <- w/Winf
    delta <- diff(c(w, Winf))
    u <- 10
    psi_F <- (1 + (w / Wfs)^-u )^-1
    psi_m <- (1 + (w / (eta_m * Winf))^-u )^-1
    m <- a * A * w^(n - 1) +  F * psi_F
    m[M] <- Inf
    g <- A*w^n *(1 - (w/Winf)^(1-n) * (epsilon_a + (1 - epsilon_a) * psi_m))
    N <- exp(- cumsum(m / g * delta)) / (g * w_r)
    N[M] <- 0
    B <- sum(psi_m  * N * w * delta)
    Rrel <- 1 - Winf^(1-n) * w_r / (epsilon_r * (1 - epsilon_a) * A * B)
    F * Rrel * sum(psi_F * N * w * delta)
  
  }
getY2 <- function(F, p = parameters())
  {
    Winf <- exp(p@logWinf) * p@scaleWinf
    A <- exp(p@logA) * p@scaleA
    n <- exp(p@logn) *p@scalen
    eta_F <- exp(p@logeta_F) * p@scaleeta_F
    eta_m <- exp(p@logeta_m) * p@scaleeta_m
    a <- exp(p@loga) *  p@scalea
    epsilon_a <- exp(p@logepsilon_a) * p@scaleepsilon_a
    epsilon_r <- exp(p@logepsilon_r) * p@scaleepsilon_r
    M <- 300
    w_r <- w_egg <- 0.001
    Delta <- (log(Winf) - log(w_r)) / (M - 1)
    w <- exp(log(w_r) + (1:M - 1) * Delta)
    wrel <- w/Winf
    delta <- diff(c(w, Winf))
    u <- 10
    psi_F <- (1 + (w / (eta_F * Winf))^-u )^-1
    psi_m <- (1 + (w / (eta_m * Winf))^-u )^-1
    m <- a * A * w^(n - 1) +  F * psi_F
    m[M] <- Inf
    g <- A*w^n *(1 - (w/Winf)^(1-n) * (epsilon_a + (1 - epsilon_a) * psi_m))
    N <- cumprod(exp(Delta)^(-m/g*w))/(g*w_r)
    N[M] <- 0
    B <- sum(psi_m  * N * w * delta)
    Rrel <- 1 - Winf^(1-n)/(epsilon_r * (1 - epsilon_a) * A * B) * (w_r/w_egg)^(a-1)
    F * Rrel * sum(psi_F * N * w * delta)
  
  }
getYR <- function(F, p = parameters())
  {
    Winf <- exp(p@logWinf) * p@scaleWinf
    A <- exp(p@logA) * p@scaleA
    n <- exp(p@logn) *p@scalen
    eta_F <- exp(p@logeta_F) * p@scaleeta_F
    eta_m <- exp(p@logeta_m) * p@scaleeta_m
    a <- exp(p@loga) *  p@scalea
    epsilon_a <- exp(p@logepsilon_a) * p@scaleepsilon_a
    epsilon_r <- exp(p@logepsilon_r) * p@scaleepsilon_r
    M <- 1000
    w_r <- w_egg <- 0.001
    Delta <- (log(Winf) - log(w_r)) / (M - 1)
    w <- exp(log(w_r) + (1:M - 1) * Delta)
    rm(Delta)
    wrel <- w/Winf
    delta <- diff(c(w, Winf))
    u <- 10
    psi_F <- (1 + (w / (eta_F * Winf))^-u )^-1
    psi_m <- (1 + (w / (eta_m * Winf))^-u )^-1
    m <- a * A * w^(n - 1) +  F * psi_F
    m[M] <- Inf
    g <- A*w^n *(1 - (w/Winf)^(1-n) * (epsilon_a + (1 - epsilon_a) * psi_m))
    N <- exp(- cumsum(m / g * delta)) / (g * w_r)
    N[M] <- 0
    B <- sum(psi_m  * N * w * delta)
    Rrel <- 1 - Winf^(1-n)/(epsilon_r * (1 - epsilon_a) * A * B) * (w_r/w_egg)^(a-1)
    F * Rrel * sum(psi_m * N * w * delta)
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

getRandomParameters.fixedWinf <- function(winf, Rrel.gt=-Inf, Fmsy.gt=0) {
  parameter.names <- c("A", "n" ,"eta_m","eta_F", "a" ,"Fm","Winf","epsilon_a", "epsilon_r")
  parameter.value <- c(4.5,0.75 , 0.25  ,  0.05 , 0.35,0.25,  winf ,    0.8    ,     0.1    )
  getRandomParameters(parameter.names, parameter.value,, Rrel.gt=Rrel.gt, Fmsy.gt=Fmsy.gt)
}

tmclapply <- function(X, FUN, ..., simplify=FALSE){
  aplfun <- if(require(multicore)) mclapply else lapply
  start <- Sys.time()
  pb <- txtProgressBar(min = 0, max = 100, style=3)
  results <- local({
    f <- fifo(tempfile(), open="w+b", blocking=TRUE)
    if (inherits(fork(), "masterProcess")) {
      progress <- 0.0
      while(progress < 1 && !isIncomplete(f)) {
        msg <- readBin(f, "double")
        progress <- progress + as.numeric(msg)
        setTxtProgressBar(pb, progress * 100)
        tt <- (1 - progress)*(difftime(Sys.time(), start, units="mins"))/ progress
        cat(" ETC:", as.integer(tt), "min(s) and", round((tt - as.integer(tt)) * 60, 0) ,"secs")
      } 
      exit()
    }
    res <- aplfun(X, function(x) {rr <- FUN(x); writeBin(1/length(X), f); rr})
    close(f)
    setTxtProgressBar(pb,100)
    close(pb)
    res
  })
  cat(difftime(Sys.time(), start, units="mins"), "mins\n")
  if (simplify) simplify2dataframe(results) else results
} 

simplify2dataframe <- function(dd) {
  data.frame(t(simplify2array(dd)))
}
