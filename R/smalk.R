weightAtLength <- function(smalkFile, years, plot=FALSE) {
  smalk <- 
  nyr <- length(years)
  res <- matrix(nrow=nyr, ncol=2, dimnames=list(years, c("a", "b")))
  if(plot) {
   setParMfrow(years)
  }
  for(y in years) {
    d <- smalk[smalk$Year == y, ]
    fitab <- lm(log(d$IndividualWeight) ~ log(d$LngtClasMM/10))
    res[as.character(y), "a"] <- exp(fitab$coefficients[1])
    res[as.character(y), "b"] <- fitab$coefficients[2]
    if(plot) {
      plot(d$LngtClasMM / 10, d$IndividualWeight)
      x <- seq(min(d$LngtClasMM / 10), max(d$LngtClasMM / 10))
      y <- exp(fitab$coefficients[1]) * x ^ fitab$coefficients[2]
      title(paste("a = ", round(exp(fitab$coefficients[1]),digits=5), ", b = ", round(fitab$coefficients[2], digits=5)))
      lines(x,y, col=2)
    }
  }
  res
}


## ## Test
## d <- readExchange("~/Work/mainCode/R/SecondPaper/cod/NorthSea/data/survey/ExchangeDataNS-IBTS1965-2014All.zip")
## cod <- subset(d, Species=="Gadus morhua")
## cod <- addSpectrum(cod)
## size <- as.data.frame(cod, format="long")

## str(size)
