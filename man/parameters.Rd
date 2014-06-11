\name{parameters}
\alias{parameters}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
  Constructor for the Parameters Class
}
\description{
  Easier constructor for the Parameters class.
}
\usage{
parameters(names = c(), vals = c())
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{names}{
    String vector. Contains the names of the parameters that will have
    non default values.
  }
  \item{vals}{
    Numeric vector. The corresponding values.
  }
}
\details{
  The values in 'vals' should be given as logarithms of the parameter values.
}
\value{
  Returns an object of the Parameters class
}

\author{
  alko
}

\seealso{
  \code{\link{Parameters-class}}, ~~~
}
\examples{
\dontrun{
## Without any arguments gives a Parameters object with default values
parameters()

## Changing some parameters gives the corresponding object
parameters(c("Winf", "Fm"), c(log(1000), log(0.4)))

}
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{constructor }