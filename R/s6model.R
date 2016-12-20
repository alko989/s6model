#' Data poor stock assessment using a single species size based model in steady state
#'
#' @references Kokkalis, A., Nielsen, A., Thygesen U. H. and K. H. Andersen. (In review) Reliability of fisheries reference points estimation for data-poor stocks, Fisheries research
#' @references Andersen, K. H., & Beyer, J. E. (2013). Size structure, not metabolic scaling rules, determines fisheries reference points. Fish and Fisheries. doi:10.1111/faf.12042
#' 
#' @docType package
#' @name s6model
#' @importFrom numDeriv jacobian hessian
#' @import methods
#' @importFrom truncdist rtrunc
#' @importFrom stats simulate acf aggregate coefficients density integrate lm na.omit nlminb optimize pchisq ppois predict qnorm qpois qqline qqnorm quantile rlnorm rmultinom rnorm runif setNames shapiro.test t.test
#' @importFrom utils packageDescription packageName head tail
#' @importFrom graphics abline axis box hist legend lines matplot mtext par plot plot.default points polygon text title
#' @importFrom grDevices col2rgb grey rgb xy.coords
#' @importFrom TMB MakeADFun newtonOption config sdreport
NULL
