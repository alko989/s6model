#' @useDynLib s6model
#' @useDynLib s6modelts
#' @useDynLib calcFmsy

checkTMBPackageVersion <- function() {
  ## Taken from kaskr/adcomp/TMB
  file <- paste0(system.file(package="s6model"),"/TMB-version")
  cur.TMB.version <- as.character(packageVersion("TMB"))
  if (!file.exists(file)) {
    writeLines(cur.Matrix.version, con = file)
  }
  s6.TMB.version <- readLines(file)
  if (!identical(s6.TMB.version, cur.Matrix.version)) {
    warning(
      "Package version inconsistency detected.\n",
      "s6model was built with TMB version ",
      s6.TMB.version,
      "\n",
      "Current Matrix version is ",
      cur.TMB.version,
      "\n",
      "Please re-install 's6model' from source using remotes::install_github('alko989/s6model').")
  }
}


.onLoad <- function(lib, pkg) {
  checkTMBPackageVersion()
  packageStartupMessage("Loading ", getVersion(),"\n")
}



.onUnload <- function(lib) {
  library.dynam.unload("s6model", lib)
  library.dynam.unload("s6modelts", lib)
  library.dynam.unload("calcFmsy", lib)
}

