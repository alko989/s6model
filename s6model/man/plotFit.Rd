\name{plotFit}
\alias{plotFit}
\alias{plotFit-methods}
\alias{plotFit,Parameters,numeric,logical-method}
\alias{plotFit,Parameters,numeric,missing-method}
\title{
  Visualisation of the fit to the data
}
\description{
  Creates a plot to visualise how good the estimated parameters fit the
  data. A kernel density estimation of the data and the pdf using
  the fitted parameters are ploted.
  %%  ~~ A concise (1-5 lines) description of what the function does. ~~
}
\usage{
plotFit(object, data, add)

}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{object}{
    A Parameters object 
  }
  \item{data}{
    Numeric vector. Contains data (i.e weights in grams of fished individuals).
  }
  \item{add}{
    Boolean. If TRUE it is added using lines to an active plot. A new
    plot is produced if FALSE (default).
  }
}

\value{
  None
}

\author{
  alko
}

\examples{
\dontrun{
plotFit(parameters(), simulateData())
}
}
\keyword{plot}

