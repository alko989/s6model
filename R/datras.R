##' Subsets a DATRAS object
##'
##' Selects species, gear, years, haul duration from DATRAS object and returns it
##' or save it in an RDS file
##' @param dat DATRAS object, or character containing the filename of an RDS file holding a DATRASraw object
##' @param species character vector, species to keep (latin name). If \code{NULL} all species are kept
##' @param gear character vector, gears to keep. If \code{NULL} all gears are kept
##' @param years numeric vector, if length is 2 it is the first and last year to keep, if length is 1 only the data of this year are kept. If \code{NULL} all years are kept
##' @param haulDur numeric vector of length 2 with the minimum and maxumum haul duration. If \code{NULL} all haul durations are kept
##' @param out string, filename of the output RDS file
##' @return if \code{out} is NULL, a DATRASraw object (invisible). If \code{out} is character, the filename of the saved RDS file. If the file specified in \code{out} exist, the function does not rewrite the file and returns the filename.
##' @author alko
##' @export
subsetDatras <- function(dat, species="Gadus morhua", gear="GOV",
                          years=c(1991, 2013), haulDur=c(25, 35), out=NULL) {
  pnh <- function(o) cat("Number of hauls: ", nrow(o[[2]]), "\n")
  if(! is.null(out)) {
    if(file.exists(out)) return(out)
  }
  if(class(dat) == "character") dat <- readRDS(dat)
  if(class(dat) != "DATRASraw") stop("Not a DATRAS object, please use the DATRAS",
                                     "package to read the exchange files")
  
  if(! species %in% levels(dat[["HL"]]$Species))
    stop("The species ", species, " is not in the provided DATRAS object")
  
  cat("\n *** Selecting the data for species: ", species, "***\n\n")
  res <- subset(dat, Species == species)
  pnh(res)
  if( ! is.null(gear)) {
    cat("\n *** Selecting the data for gear: ", gear, "***\n\n")
    res <- subset(res, Gear %in% gear)
    pnh(res)
  }
  if( ! is.null(years)) {
    cat("\n *** Selecting the data from the years: ", paste(years, collapse=" - "), "***\n\n")
    if(length(years) == 1)
      res <- subset(res, Year == years)
    if(length(years) == 2)
      res <- subset(res, Year %in% seq(years[1], years[2]))
    if(length(years) == 0 | length(years) >2) stop("years has to be a numeric vector of length 1 or 2")
    pnh(res)
  }
  if( ! is.null(haulDur)) {
    cat("\n *** Selecting only hauls with duration: ", paste(haulDur, collapse=" - "), "***\n\n")
    res <- subset(res, haulDur[1]<HaulDur & HaulDur<haulDur[2])
    pnh(res)
  }
  cat("\n *** Adding spectrum ***\n\n")
  res <- addSpectrum(res)
  if(! is.null(out)) {
    saveRDS(res, file=out)
    cat("\nSaved data in: ", out, "\n\n")
    return(out)
  }  
  invisible(res)   
}

##' Get yearly weight frequencies
##'
##' For selected years get yearly weight frequencies from DATRAS object
##' @param dat DATRASraw object
##' @param years numeric vector, containing the years to return
##' @param ... additional arguments to \code{\link{datrasraw2weightfreq}}
##' @return list of \code{data.frame}s, one for each year
##' @author alko
##' @export
getDfYears <- function(dat, years = as.numeric(levels(dat[[2]]$Year)), ...) {
  lapply(years, function(yr) {
    d <- subset(dat, Year %in% yr)
    datrasraw2weightfreq(d, ...)
  })
}

##' Convert a DATRASraw object to weight frequency data.frame
##'
##' The length information of a DATRASraw object is transformed to weight and the weight frequency is returned as \code{data.frame}
##' @param datr DATRASraw object
##' @param a numeric weight-length relationship, parameter a, see Details
##' @param b numeric, weight-length relationship, parameter b, see Details
##' @param estWL logical, if TRUE and available information about weight is available, the weight-length relationship is fitted
##' @param verbose logical, if TRUE it shows information about the fitted weight length parameters and the output \code{data.frame}
##' @return data.frame with columns Weight, Freq containing the weight frequencies of the input data
##' @details The weight-length relationship has the form $W = aL^b$
##' @author alko
##' @export
datrasraw2weightfreq <- function(datr, a=0.01, b=3, estWL=FALSE, verbose=TRUE) {
  
  df <- aggregate(Count  ~ LngtCm, data=datr[["HL"]], sum, na.action=na.omit)
  names(df) <- c("Length", "Freq")
  if(estWL) {
    fit <- lm(log(IndWgt) ~ log(LngtCm), data = datr[[1]])
    coef <- coefficients(fit)
    if(verbose) cat(coef, "\n")
    a <- exp(coef[1])
    b <- coef[2]
    if(verbose) print(summary(fit))
  }
  df$Weight <- a * df$Length ^ b
  if(verbose) showDf(df)
  attr(df, "createdBy") <- getVersion()
  invisible(df)
}

showDf <- function(df) {
  if(dim(df)[1] > 12) {
    print(head(df,4), digits = 2)
    cat("+------------------------+\n",
        "   Ommited ", dim(df)[1] - 8, " lines\n",
        "+------------------------+\n", sep="")
    print(tail(df, 4),digits = 2)
  } else {
    print(df, 2)
  }
}




