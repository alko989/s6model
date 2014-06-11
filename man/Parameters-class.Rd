\name{Parameters-class}
\Rdversion{1.1}
\docType{class}
\alias{Parameters-class}
\alias{lines,Parameters-method}
\alias{plot,Parameters-method}
\alias{show,Parameters-method}

\title{Class \code{"Parameters"}}
\description{
 A class that contains all parameters of the single species size-based model.
}
\section{Objects from the Class}{
Objects can be created by calls of the form \code{new("Parameters", ...)}.
Alternatively the 'constructor' \code{\link{parameters}} can be used.
}
\section{Slots}{
  \describe{
    \item{\code{logWinf}:}{Object of class \code{"numeric"} ~~ }
    \item{\code{logFm}:}{Object of class \code{"numeric"} ~~ }
    \item{\code{logh}:}{Object of class \code{"numeric"} ~~ }
    \item{\code{logn}:}{Object of class \code{"numeric"} ~~ }
    \item{\code{logeta_F}:}{Object of class \code{"numeric"} ~~ }
    \item{\code{logeta_m}:}{Object of class \code{"numeric"} ~~ }
    \item{\code{loga}:}{Object of class \code{"numeric"} ~~ }
    \item{\code{logepsilon_a}:}{Object of class \code{"numeric"} ~~ }
    \item{\code{logepsilon_r}:}{Object of class \code{"numeric"} ~~ }
    \item{\code{M}:}{Object of class \code{"numeric"} ~~ }
  }
}
\section{Methods}{
  \describe{
    \item{lines}{\code{signature(x = "Parameters")}: ... }
    \item{plot}{\code{signature(x = "Parameters")}: ... }
    \item{show}{\code{signature(object = "Parameters")}: ... }
	 }
}

\author{
alko
}


\seealso{
\code{\link{parameters}}

}
\examples{
showClass("Parameters")
}
\keyword{classes}
