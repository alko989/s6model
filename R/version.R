getVersion <- function(pkg=packageName()) {
  pd <- packageDescription(pkg)
  v <- paste0(pd$Package, "_v", pd$Version)
  if(is.null(pd$GithubRef))
    return(v)
  else
    paste0(v , "@", pd$GithubSHA1)
}

#' Adds version watermark to a plot
#' 
#' Adds the version and, if installed from Github, the first 10 characters from
#' SHA commit code
#' 
#' @param cex numeric, controls the font size
#' @param col A color code or name
#' @return Invisible NULL
#' @seealso "Color specification" in the documentation of \code{\link{par}}
#' @author alko
#' @keywords aplot
#' @examples
#' 
#' hist(simulateData3()$sample)
#' addVersion()
#' 
#' @export
addVersion <- function(cex=0.5, col="#12345655", lengthSHA = 6) {
    v <- getVersion()
    mtext(formatVersion(v, lengthSHA = lengthSHA), side=4, line=-0.3, adj=0.01, col = col, cex = cex)
}

formatVersion <- function(x, showSHA = TRUE, lengthSHA = 6) {
  pkgname <- regmatches(x, regexpr(pattern = "^[^_]*", v, perl = TRUE))
  ver <- regmatches(x, regexpr(pattern = "v[0-9.]*", v, perl = TRUE))
  res <- paste(pkgname, ver)
  if(grepl("@", x) & showSHA) {
    sha <- regmatches(x, regexpr(pattern = paste0("[^@]{",lengthSHA, "}$"), v, perl = TRUE))
    res <- paste0(res, " (@", sha, ")")
  }
  res
  
}
