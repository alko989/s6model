#' @useDynLib s6model
.onLoad <- function(lib, pkg) {
  library.dynam("s6model", pkg, lib)
}

.onUnload <- function (lib) {
  library.dynam.unload("s6model", lib)
}
