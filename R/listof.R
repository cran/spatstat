#
# listof.R
#
# Methods for class `listof'
#
# plot.listof is defined in plot.splitppp.R
#

"[<-.listof" <- function(x, i, value) {
  # invoke list method
  class(x) <- "list"
  x[i] <- value
  # then make it a 'listof' object too
  class(x) <- c("listof", class(x))
  x
}
  
summary.listof <- function(object, ...) {
  x <- lapply(object, summary, ...)
  class(x) <- "summary.listof"
  x
}

print.summary.listof <- function(x, ...) {
  class(x) <- "listof"
  print(x)
  invisible(NULL)
}

listof <- function(...) {
  stuff <- list(...)
  class(stuff) <- c("listof", class(stuff))
  return(stuff)
}

as.listof <- function(x) {
  if(!is.list(x))
    x <- list(x)
  class(x) <- c("listof", class(x))
  return(x)
}

contour.listof <- function(x, ...) {
  plot(x, ..., plotcommand="contour")
}

image.listof <- function(x, ...) {
  plot(x, ..., plotcommand="image")
}
