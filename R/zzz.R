#' @useDynLib s6model
#' @useDynLib s6modelts
#' @useDynLib calcFmsy
.onAttach <- function(lib, pkg) {
   # packageStartupMessage("Loading ", s6model:::getVersion(),"\n")  
 }

.onUnload <- function (lib) {
  library.dynam.unload("s6model", lib)
  library.dynam.unload("s6modelts", lib)
  library.dynam.unload("calcFmsy", lib)
}

