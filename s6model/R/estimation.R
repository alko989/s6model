reloadSteadyStateEstimation <- function() {
  detach("package:steadyStateEstimation", unload=TRUE)
  library(steadyStateEstimation)
}
    

minimizeme <- function(theta, data, names, fixed.names=c(), fixed.vals=c(), isSurvey=FALSE)
{
  params <- parameters(c(names, fixed.names), c(theta, fixed.vals))  
  if(class(data) == "data.frame") {
    
    return(with(getParams(params,isSurvey),
                sum( - data$Freq * log(pdfN.approx(as.numeric(as.character(data$Weight)) )))))
  }
  return(with(getParams(params,isSurvey),
              sum(-log(pdfN.approx(data)))))
}

estimateParam <-
  function(names = c("Fm", "a"),
           data=simulateData3(parameters(), samplesize=1000)$sample,
           bootstrap=1, start= rep(0.5, length(names)),
           lower = rep(-Inf, length(names)), upper=rep(Inf, length(names)),
           fixed.names=c(), fixed.vals=numeric(0),
           plotFit=FALSE, isSurvey=FALSE, verbose=getOption("verbose"), ...)
  {
    p <- parameters()
    if(is.null(lower))
        lower <- rep(-Inf, length(names))
    lower[which(names == "Winf")] <- 
        if(class(data)=="data.frame") 
            log((max(data$Weight) + 1) / p@scaleWinf)
        else
            log((max(data) + 1) / p@scaleWinf)
        
    scales <- sapply(names, function(n) get(paste0("getscale", n))(p))
    
    if(require(multicore))
    {
      useapply <- mclapply
    } else {
      useapply <- lapply
    }
    
    sd <- mean <- rep(0, length(names))
    if(bootstrap == 1)
      {
        estim <- nlminb(log(start), minimizeme, data=data, names=names,
                        fixed.names=fixed.names, fixed.vals=fixed.vals,
                        lower=lower, upper=upper, isSurvey=isSurvey)
        
        if(estim$convergence != 0) {
          warning(estim$message)
        }
        res <- c(estim$par) ##, estim$convergence)
        h <- numDeriv::hessian(minimizeme, estim$par, data=data,
                     names=names, fixed.names=fixed.names,
                     fixed.vals=fixed.vals,isSurvey=isSurvey)
        s <- numDeriv::jacobian(minimizeme, estim$par, data=data,
                     names=names, fixed.names=fixed.names,
                     fixed.vals=fixed.vals,isSurvey=isSurvey)
        vcm <- try(solve(h))
        if(class(vcm) != "try-error") {
          st.er <- sqrt(diag(vcm)) 
          ci <- cbind(exp(res)*scales,
                      exp(outer(1.96 * st.er, c(-1,1), '*') + res) * scales)
        } else {
          st.er <- NA
          ci <- matrix(rep(NA, length(names)*3),ncol=3)
        }
        
        dimnames(ci) <- list(names, c("Estimate","Lower", "Upper"))
      }
    else
      {
        res <- simplify2array(useapply(1:bootstrap, function(x) {
          estim <-
            nlminb(log(start), minimizeme, data=sample(data, replace=TRUE),
                   names=names, lower=log(lower), fixed.names=fixed.names,
                   fixed.vals=fixed.vals, isSurvey=isSurvey, Frasier=Frasier)
          if(estim$convergence != 0)
            warning(estim$message)
          estim$par##, convergence = estim$convergence)
        }))
        if(length(names) == 1){
          rr <- exp(res)*scales
          sd <- sd(rr)
          mean <- mean(rr)
          res <- log(mean / scales)
        } else {
          rr <- sapply(1:bootstrap, function(y) exp(res[ ,y]) * scales)
          sd <- apply(rr, 1, sd)
          mean <- apply(rr, 1, mean)
          res <- log(mean / scales)
        }
      }
    
    p <- parameters(c(names, fixed.names), c(t(simplify2array(res)), fixed.vals))

    if(plotFit) { plotFit(p, data, ...) }
    if(bootstrap == 1) {
      if(verbose)
          print(ci)
      return(structure(p, par=res, hessian=h, jacobian=s, st.er=st.er, ci=ci, objective=estim$objective, convergence=estim$convergence, nlminbMessage=estim$message, call=match.call()))
    } else {
      return(list(parameters = p, mean=mean, sd=sd,allResults=rr, names=names))
    }
  }

minimizemeMultidata <- function(theta, surdata, comdata, names, fixed.names=c(), fixed.vals=c(),  Frasier=FALSE)
  {
    params <- parameters(c(names, fixed.names), c(theta, fixed.vals))
    return(with(getParams(params,isSurvey=TRUE, surveySelectivityFrasier=Frasier), sum(-log(pdfN.approx(surdata)))) +
           with(getParams(params,isSurvey=FALSE, surveySelectivityFrasier=Frasier), sum(-log(pdfN.approx(comdata)))))
  }

estimateMultidata <-
  function(names = c("Fm", "Winf", "eta_F"),
           surdata=simulateData3(parameters(), samplesize=1000, isSurvey=TRUE)$sample,
           comdata=simulateData3(parameters(), samplesize=1000)$sample,
           bootstrap=1, start= rep(0.5, length(names)),
           lower = rep(0.1, length(names)),
           fixed.names=c(), fixed.vals=numeric(0),
           plotFit=FALSE, Frasier=FALSE,...)
  {
    lower[which(names == "Winf")] <- (max(surdata, comdata) + 1) / parameters()@scaleWinf
    lower[which(names == "eta_F")] <- 0.0001
    if(require(multicore))
    {
      useapply <- mclapply
    } else {
      useapply <- lapply
    }
    sd <- mean <- seq(along.with=names)*0
    if(bootstrap == 1)
      {
        estim <- nlminb(log(start), minimizemeMultidata, surdata=surdata, comdata=comdata,
                        names=names, fixed.names=fixed.names, fixed.vals=fixed.vals,
                        lower=log(lower), Frasier=Frasier)
        if(estim$convergence != 0)
          warning(estim$message)
        res <- c(estim$par) ##, estim$convergence)
      }
    else
      {
        res <- simplify2array(useapply(1:bootstrap, function(x) {
          estim <-
            nlminb(log(start), minimizeme, data=sample(data, replace=TRUE),
                   names=names, lower=log(lower), fixed.names=fixed.names,
                   fixed.vals=fixed.vals, isSurvey=isSurvey, Frasier=Frasier)
          if(estim$convergence != 0)
            warning(estim$message)
          estim$par##, convergence = estim$convergence)
        }))
        if(length(names) == 1){
          scale <- eval(parse(text=paste("parameters()@scale",names,sep="")))
          rr <- exp(res)*scale
          sd <- sd(rr)
          mean <- mean(rr)
          res <- log(mean / scale)
        } else {
          scales <- sapply(seq(along.with=names), function(x)
                           eval(parse(text=paste("parameters()@scale",
                                        names[x],sep=""))))
          rr <- sapply(1:bootstrap, function(y) exp(res[ ,y]) * scales)
          sd <- apply(rr, 1, sd)
          mean <- apply(rr, 1, mean)
          res <- log(mean / scales)
        }
      }
    
    p <- parameters(c(names, fixed.names), c(t(simplify2array(res)), fixed.vals))
    if(plotFit) plotFit(p, data, ...)
    if(bootstrap == 1)
      return(p)
    else
      return(list(parameters = p, mean=mean, sd=sd,allResults=rr, names=names))
  }



minimizeme2 <- function(theta, data, names, fixed.names=c(), fixed.vals=c())
  {
    params <- parameters(c(names, fixed.names), c(theta, fixed.vals))
    return(with(getParams2(params), sum(-log(pdfN.approx(data)))))
  }

estimateParam2 <- function(names = c("Fm", "a"),
                          data=simulateData3(parameters(), samplesize=1000)$sample,
                          bootstrap=1, start= rep(0.5, length(names)), lower = rep(0.1, length(names)),
                           fixed.names=c(), fixed.vals=numeric(0),
                          plotFit=FALSE)
  {
    lower[which(names == "Winf")] <- (max(data) + 1) / parameters()@scaleWinf
      
    if(require(multicore))
    {
      useapply <- mclapply
    } else {
      useapply <- lapply
    }
    if(bootstrap == 1)
      {
        estim <- nlminb(log(start), minimizeme2, data=data, names=names,
                        fixed.names=fixed.names, fixed.vals=fixed.vals, lower=log(lower))
        if(estim$convergence != 0)
          warning(estim$message)
        res <- c(estim$par) ##, estim$convergence)
      }
    else
      {
        res <- useapply(1:bootstrap, function(x) {
          estim <-
            nlminb(log(start), minimizeme, data=sample(data, replace=TRUE),
                   names=names, lower=log(lower), fixed.names=fixed.names, fixed.vals=log(fixed.vals))
          if(estim$convergence != 0)
            warning(estim$message)
          estim$par##, convergence = estim$convergence)
        })
      }
    ## estim <- nlminb(start, minimizeme, data=data, names=names)
    ## t(simplify2array(res))
    p <- parameters(c(names, fixed.names), c(t(simplify2array(res)), fixed.vals))
    if(plotFit) plotFit(p, data)
    return(p)
  }
