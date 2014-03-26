\name{simulateData}
\alias{simulateData, simulateData3}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
  Simulation of fishing data
}
\description{
  Given a Parameters object a sample of 'samplesize' is simulated.
}
\usage{
simulateData(samplesize = 1000, params = parameters())
simulateData3(samplesize = 1000, params = parameters(), ...)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{samplesize}{
    Integer. The sample size.
  }
  \item{params}{
    Parameters object.
  }
  \item{...}{
    Passes named arguments to getParams. Mainly to be used to simulate
    survey data using isSurvey=TRUE.
  }
}
\details{
 The parameters can be set using the parameters constructor function. 
}
\value{
  Returns a list containing: 'sample' - numeric vector of length 'samplesize' with
  simulated data, 'parameters' - Parameters object with the model parameters.
}

\author{
  alko@aqua.dtu.dk
}

\seealso{
  \code{\link{parameters}}
}
\examples{
## Generate a sample of length 1000 with default parameters.
sam<-simulateData3()$sample

\dontrun{
## Generate a larger sample (10000) with non-default parameters.
sam2<-simulateData3(samplesize=10000,
                    params=parameters(c("a", "Fm"), c(log(0.2), log(0.3))))
}
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{misc}

