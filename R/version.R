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
#' Adds the version and, if installed from Github, the first part of the
#' SHA commit code
#' 
#' @param v character vector, version(s) to add to plot
#' @param cex numeric, controls the font size
#' @param col A color code or name
#' @param lengthSHA integer, amount of characters to be printed
#' @param description character vector, same length as \code{x}, description of each version
#' @return Invisible NULL
#' @details All versions in the \code{v} vector are added to the current plot,
#' each with a discription from \code{description} if it is not \code{NULL}.
#' 
#' If \code{v} is \code{NULL}, the version returned by the internal function \code{getVersion} is used.
#' @seealso "Color specification" in the documentation of \code{\link{par}}
#' @author alko
#' @keywords aplot
#' @examples
#' 
#' hist(simulate(s6params()))
#' addVersion()
#' 
#' @export
addVersion <- function(v = NULL, cex=0.5, col="#12345655", lengthSHA = 6, description = NULL) {
  if(is.null(v)) v <- getVersion()
  fv <- sapply(v, formatVersion, lengthSHA = lengthSHA, USE.NAMES = FALSE)
  text <- if(length(description) == length(v)) {
    paste0(mapply(function(v, d) {
      paste0(v, " (",d,")")
    }, fv, description), collapse = " - ")
  } else {
    paste(fv, collapse = " - ")
  } 
  
  mtext(text, side=4, line= 0.1, padj = 0.1, adj=0.01, col = col, cex = cex)
}

formatVersion <- function(x, showSHA = TRUE, lengthSHA = 6) {
  pkgname <- regmatches(x, regexpr(pattern = "^[^_]*", x, perl = TRUE))
  ver <- regmatches(x, regexpr(pattern = "v[0-9.]*", x, perl = TRUE))
  res <- paste(pkgname, ver)
  if(grepl("@", x) & showSHA) {
    sha <- regmatches(x, regexpr(pattern = paste0("[^@]{",lengthSHA, "}$"), x, perl = TRUE))
    res <- paste0(res, " (@", sha, ")")
  }
  res
  
}
