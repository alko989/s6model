\name{getParams}
\alias{getParams}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
 Returns a list with all parameters and functions
}
\description{
  Given a Parameters object, it returns as a list all relevant parameters and
  corresponding functions like pdf and cdf.  
}
\usage{
getParams(p = new("Parameters"), isSurvey=FALSE, surveySelectivityFrasier=FALSE, calcFmsy=FALSE)
}
\arguments{
  \item{p}{
    A Parameters object
  }
  \item{isSurvey}{
    boolean, if TRUE a survey selectivity is used for the pdf
  }
  \item{surveySelectivityFraser}{
    boolean, if TRUE the selectivity described in Fraser et al(2007) for
    North Sea cod
  }
  \item{calcFmsy}{
    boolean, if TRUE the Fmsy is caculated using optimization
  }
}

\value{
  An invisible list, containing all model parameters and the cdf, pdf functions. Also
  calculations of the Spawning Stock Biomass and the Yield and if
  'calcFmsy' is TRUE, the Fmsy.
}

\author{
 alko@aqua.dtu.dk
}

\seealso{
\code{\link{parameters}}
}
\examples{
 getParams()
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{misc}
