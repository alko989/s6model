#' @useDynLib calcFmsy
#' @useDynLib s6model
use.onLoad <- function(lib, pkg) {
  library.dynam("calcFmsy", pkg, lib)
  library.dynam("s6model", pkg, lib)
}