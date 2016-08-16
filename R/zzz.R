#' @useDynLib s6model
#' @useDynLib s6modelts
#' @useDynLib calcFmsy
#' @useDynLib s6modelcomb
.onAttach <- function(lib, pkg) {
   packageStartupMessage("Loading ", getVersion(),"\n")  
 }

.onUnload <- function (lib) {
  library.dynam.unload("s6model", lib)
  library.dynam.unload("s6modelts", lib)
  library.dynam.unload("calcFmsy", lib)
  library.dynam.unload("s6modelcomb", lib)
}

