#' Data-limited stock assessment using the single-species, size-structured, steady-state (s6) model
#'
#' @references Kokkalis, A., Thygesen, U. H., Nielsen, A. & Andersen, K. H. (2015). Limits to the reliability of size-based fishing status estimation for data-poor stocks. Fisheries Research, 171, 4–11. doi:10.1016/j.fishres.2014.10.007 \href{https://dx.doi.org/10.1016/j.fishres.2014.10.007}{↗️}
#' @references Kokkalis, A., Eikeset, A. M., Thygesen, U. H., Steingrund, P. & Andersen, K. H. (2016). Estimating uncertainty of data limited stock assessments. ICES Journal of Marine Science: Journal Du Conseil, 74(1), 69–77. doi:10.1093/icesjms/fsw145 \href{https://dx.doi.org/10.1093/icesjms/fsw145}{↗}
#' @references Andersen, K. H. & Beyer, J. E. (2015). Size structure, not metabolic scaling rules, determines fisheries reference points. Fish and Fisheries, 16(1), 1-22 . doi:10.1111/faf.12042 \href{https://dx.doi.org/10.1111/faf.12042}{↗}
#'
#' @docType package
#' @name s6model
#' @importFrom numDeriv jacobian hessian
#' @import methods
#' @importFrom truncdist rtrunc
#' @importFrom stats simulate
#' @importFrom utils packageDescription
#' @importFrom utils packageName
NULL
