#' Make a report about the (quasi-)residuals
#'
#' @param rr s6output object, as returned by makeAssessment
#'
#' @return list with residuals and the results of the t-test, acf and Shapiro-Wilks test
#' @export
plotResiduals <- function(rr) {
  par(mfrow = c(3,2), mar = c(4.1,4.1, 1,1))
  obj <- attr(rr, "obj")
  rprt <- attr(rr, "obj")$env$report(par = attr(rr, "opt")$par)
  obs <- rprt$freq
  exp <- rprt$Nvec / sum(rprt$Nvec) * rr$sigma
  w <- which(obs > 0)
  # Observed count and weight
  plot(rprt$Weight[w], obs[w], pch = 20, ylim = range(obs[w], exp[w]), log = "xy", col = 'lightgrey',xlab = 'Weight (g)', ylab = 'Observed count') # sub = 'Size distribution', 
  lines(rprt$Weight[w], exp[w])
  D <- 2 * ( obs[w] * log(obs[w] / exp[w]) - (obs[w] - exp[w]))
  chi.sq <- (obs[w]-exp[w])^2/exp[w]
  degf <- length(w) - 3 - 1
  1 - pchisq(sum(chi.sq), degf)
  1 - pchisq(sum(D), degf)
  
  plot(exp[w], obs[w],pch = 20, cex = 0.7, ylab = "Observed count", xlab = "Expected count") #, sub = 'Observed vs. expected count'
  mx <- max(obs[w], exp[w]) * 1.1
  lower.quantile <- function(x)qpois(.025,x)
  plot(lower.quantile,0,mx,add=TRUE,col="lightgrey",n=1e4)
  upper.quantile <- function(x)qpois(.975,x)
  plot(upper.quantile,0,mx,add=TRUE,col="lightgrey",n=1e4)
  abline(0,1, col = "lightgrey")
  perc.in <- round(sum(exp[w] >= lower.quantile(obs[w]) & exp[w] <= upper.quantile(obs[w])) / length(w) * 100, 2)
  
  # Pseudo residuals for poisson (see Zuchinni book)
  rsd <- ppois(obs[w], lambda = exp[w])
  qrsd <- qnorm(rsd)
  qrsd <- qrsd[is.finite(qrsd)]
  if (length(qrsd) < length(rsd)) cat("Some of the residuals are not finite (", length(rsd) - length(qrsd)," out of ", length(rsd), ")", sep ="")
  plot(qrsd, ylab = 'Pseudo-residuals') #  sub = 'Residuals',
  qqnorm(qrsd, main = NULL)
  qqline(qrsd, col= 'lightgrey')
  
  ## Normality test (Ho: normaly distributed, p > 0.05 => fail to reject)
  st <- shapiro.test(qrsd)
  isNormal <- st$p.value > 0.05
  ## Bias test (Ho: mean of pseudo-residuals is equal to 0)
  tt <- t.test(qrsd)
  isBiased <- tt$p.value < 0.05
  ## Check for correlations
  ac <- acf(qrsd) #, main = 'Autocorrelation'
  ci.acf <- qnorm((1 + 0.95) / 2) / sqrt(length(qrsd))
  sigAc <- any(unlist(ac$acf)[-1] > ci.acf | unlist(ac$acf)[-1] < -ci.acf)
  
  msg <- "The pseudo-residuals "
  msgsp <- paste0(rep(" ", nchar(msg)), collapse = "")
  txt <- paste0("Report about residuals: \n", msg, "are ",  if(!isNormal) "not ", "normal (P = " , round(st$p.value, 4), ")\n",
               msgsp, "are ", if(!isBiased) "not ", "biased (P = " , round(tt$p.value, 2), ")\n",
               msgsp, if(!sigAc) "do not ", "have significant lags in ACF.\n\n",
               perc.in, "% of the obesrved counts are inside the\n\texpected Poisson bands\n", collapse = NULL)
  par(mar = c(0,0,0,0))
  plot(1,type = "n", axes = FALSE, xlim = c(0,1), ylim = c(0,1), xaxs = "i", yaxs = "i")
  text(0.05, 0.6, labels = txt, cex = 1, adj = 0)
  return(list(acf=ac, ttest = tt, swtest = st, residuals = qrsd))
}