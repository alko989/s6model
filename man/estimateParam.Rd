\name{estimateParam}
\alias{estimateParam}
\title{Estimates parameters using the given data}
\usage{
  estimateParam(names = c("Fm", "a"),
    data = simulateData(parameters(c("Fm","a"), log(c(0.2,
0.4))),samplesize = 1000), 
    bootstrap = 1, start = rep(log(0.5), length(names)),
    lower = log(0.01), fixed.names = c(), fixed.vals = c(),
    plotFit = TRUE)
}
\arguments{
  \item{names}{String value. The parameters to be
  estimated}

  \item{data}{Numeric string. Real or simulated catch
  length data}

  \item{bootstrap}{Integer. If it is equal 1, makes the
  estimation using the original data. If it is greater than
  one resamples the data with replacement and estimates the
  parameters 'bootstrap' times.}

  \item{start}{Numeric vector. Initial values of the
  parameters.}

  \item{lower}{The lower bound for the parameter
  estimation.}

  \item{fixed.names,fixed.vals}{String vector and numeric
  vector. Contain the names and correspondig values of
  parameters not to be estimated.}

  \item{plotFit}{Boolean. If TRUE a plot with the fited pdf and the density
    kernel estimator is produced.}
}
\value{
  The estimated values for the parameters if 'bootstrap' ==
  1. The estimated values for each resampled dataset if
  'bootstrap' is > 1.
}
\description{
  Parameter estimation minimizing the negative log
  likelihood, using the nlminb function
}
\author{
  alko
}


\examples{
## Simulate some data
sam <- simulateData(params=parameters("a", log(0.5)))

## Estimate the a parameter and see the fitted plot
estimateParam(names="a", data=sam)
}

\keyword{optimize}
