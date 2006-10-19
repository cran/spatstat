#
# Functions for extracting and setting the name of the unit of length
#
#   $Revision: 1.7 $   $Date: 2006/10/18 05:25:15 $
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


  
as.units <- function(s) {
  if(is.null(s))
    return(c("unit", "units"))
  if(!is.character(s))
    stop("unit name should be a character string or strings")
  if(length(s) == 1)
    s <- rep(s, 2)
  else if(length(s) != 2)
    stop("unit name should be a single character string, or 2 strings")
  return(s)
}
