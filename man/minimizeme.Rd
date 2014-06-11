\name{minimizeme}
\alias{minimizeme}
\title{Calculates the negative log likelihood for given parameters}
\usage{
  minimizeme(theta, data, names, fixed.names = c(),
    fixed.vals = c())
}
\arguments{
  \item{theta}{Numeric vector of parameter values}

  \item{data}{Numeric vector. The observed values, fished
  individuals in grams}

  \item{names}{String vector. Contains the names of the
  estimated parameters.}

  \item{fixed.names,fixed.vals}{String vector and numeric
  vector. Contain the names and correspondig values of
  parameters not to be estimated.}
}
\value{
  The negative log likelihood for the given parameters
}
\description{
  Is used by a minimizing function to find the parameter
  values that minimize negative log likelihood.
}
\author{
  alko
}

\examples{
\dontrun{
## We simulate some data with all default parameter values and Fm=0.4
sim <- simulateData(params=parameters("Fm", log(0.4)))

## We can plot the negative log likelihood for different values of Fm
s <- seq(0.1, 1, 0.001)
plot(sapply(s, function(x) minimizeme(theta = x, data = sim, names = "Fm"))) 

## The function is currently defined as
function (theta, data, names, fixed.names = c(), fixed.vals = c()) 
{
    params <- parameters(c(names, fixed.names), c(theta, fixed.vals))
    return(with(getParams(params), sum(-log(pdfN.approx(data)))))
}
}
}

\keyword{optimize }

