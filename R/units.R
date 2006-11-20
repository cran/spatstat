#
# Functions for extracting and setting the name of the unit of length
#
#   $Revision: 1.11 $   $Date: 2006/11/20 04:58:54 $
#
#

units <- function(x) {
  UseMethod("units")
}

units.owin <- function(x) {
  u <- as.units(x$units)
  return(u)
}

units.ppp <- function(x) {
  u <- as.units(x$window$units)
  return(u)
}

units.im <- function(x) {
  u <- as.units(x$units)
  return(u)
}

units.default <- function(x) {
  return(as.units(attr(x, "units")))
}

"units<-" <- function(x, value) {
  UseMethod("units<-")
}

"units<-.owin" <- function(x, value) {
  x$units <- as.units(value)
  return(x)
}

"units<-.ppp" <- function(x, value) {
  w <- x$window
  units(w) <- value
  x$window <- w
  return(x)
}

"units<-.im" <- function(x, value) {
  x$units <- as.units(value)
  return(x)
}

"units<-.default" <- function(x, value) {
  attr(x, "units") <- as.units(value)
  return(x)
}


###  class 'units'

makeunits <- function(sing="unit", plur="units", mul = 1) {
  if(!is.character(sing))
    stop("First entry should be a character string")
  if(!is.character(plur))
    stop("Second entry should be a character string")
  if(!is.numeric(mul) || length(mul) != 1 || mul <= 0)
    stop("Third entry should be a positive number")
  u <- list(singular=sing, plural=plur, multiplier=mul)
  if(mul != 1 && (sing=="unit" || plur=="units"))
    stop(paste("A multiplier is not allowed",
               "if the unit does not have a specific name"))
  class(u) <- "units"
  return(u)
}
  
as.units <- function(s) {
  s <- as.list(s)
  n <- length(s)
  if(n > 3)
    stop(paste("Unit name should be a character string,",
               "or a vector/list of 2 character strings,",
               "or a list(character, character, numeric)"))
  
  out <- switch(n+1,
                makeunits(),
                makeunits(s[[1]], s[[1]]),
                makeunits(s[[1]], s[[2]]),
                makeunits(s[[1]], s[[2]], s[[3]]))
  return(out)
}

print.units <- function(x, ...) {
  mul <- x$multiplier
  if(mul == 1)
    cat(paste(x$singular, "/", x$plural, "\n"))
  else 
    cat(paste(mul, x$plural, "\n"))
  return(invisible(NULL))
}
            
summary.units <- function(object, ...) {
  x <- object
  scaled <- (x$multiplier != 1)
  named  <- (x$singular != "unit")
  vanilla <- !named && !scaled
  out <-
    if(vanilla) {
      list(legend = NULL,
           axis   = NULL, 
           explain = NULL,
           singular = "unit",
           plural   = "units")
    } else if(named & !scaled) {
      list(legend = paste("Unit of length: 1", x$singular),
           axis   = paste("(", x$plural, ")", sep=""),
           explain = NULL,
           singular = x$singular,
           plural   = x$plural)
    } else {
      expanded <- paste(x$multiplier, x$plural)
      list(legend = paste("Unit of length:", expanded),
           axis   = paste("(one unit = ", expanded, ")", sep=""),
           explain  = paste("(one unit = ", expanded, ")", sep=""),
           singular = "unit",
           plural   = "units")
    }
  out <- append(out, list(scaled  = scaled,
                          named   = named,
                          vanilla = vanilla))
  class(out) <- "summary.units"
  return(out)
}

print.summary.units <- function(x, ...) {
  if(x$vanilla)
    cat("Unit of length (unnamed)\n")
  else
    cat(paste(x$legend, "\n"))
  invisible(NULL)
}



