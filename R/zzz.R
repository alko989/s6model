#' @useDynLib s6model
#' @useDynLib calcFmsy
.onAttach <- function(lib, pkg) {
   packageStartupMessage("Loading ", getVersion(),"\n")  
 }

.onUnload <- function (lib) {
  library.dynam.unload("s6model", lib)
  library.dynam.unload("calcFmsy", lib)
}

