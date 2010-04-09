#
#   distfun.R
#
#   distance function (returns a function of x,y)
#
#   $Revision: 1.3 $   $Date: 2010/04/07 08:04:59 $
#

distfun <- function(X, ...) {
  UseMethod("distfun")
}

distfun.ppp <- function(X, ...) {
  # this line forces X to be bound
  stopifnot(is.ppp(X))
  g <- function(x,y) {
    Y <- if(!missing(y)) ppp(x, y, window=X$window) else as.ppp(x)
    nncross(Y, X)$dist
  }
  attr(g, "Xclass") <- "ppp"
  class(g) <- c("distfun", class(g))
  return(g)
}

distfun.psp <- function(X, ...) {
  # this line forces X to be bound
  stopifnot(is.psp(X))
  g <- function(x,y) {
    Y <- if(!missing(y)) ppp(x, y, window=X$window) else as.ppp(x)
    nncross(Y, X)$dist
  }
  attr(g, "Xclass") <- "psp"
  class(g) <- c("distfun", class(g))
  return(g)
}

distfun.owin <- function(X, ..., invert=FALSE) {
  # this line forces X to be bound
  stopifnot(is.owin(X))
  P <- as.psp(as.polygonal(X))
  R <- as.rectangle(X)
  g <- function(x,y) {
    if(missing(y)) {
      Y <- as.ppp(x)
      x <- Y$x
      y <- Y$y
    }
    inside <- inside.owin(x, y, X)
    D <- nncross(list(x=x,y=y), P)$dist
    zero <- if(!invert) inside else !inside
    out <- ifelse(zero, 0, D)
    return(out)
  }
  attr(g, "Xclass") <- "owin"
  class(g) <- c("distfun", class(g))
  return(g)
}

print.distfun <- function(x, ...) {
  xtype <- attr(x, "Xclass")
  typestring <- switch(xtype,
                       ppp="point pattern",
                       psp="line segment pattern",
                       owin="window",
                       "unrecognised object")
  cat(paste("Distance function for", typestring, "\n"))
  X <- get("X", envir=environment(x))
  print(X)
  return(invisible(NULL))
}

plot.distfun <- function(x, ...) {
  xname <- deparse(substitute(x))
  X <- get("X", envir=environment(x))
  Z <- as.im(x, as.rectangle(X))
  do.call("plot",
          resolve.defaults(list(Z),
                           list(...),
                           list(main=xname)))
}

as.owin.distfun <- function(W, ..., fatal=TRUE) {
  X <- get("X", envir=environment(W))
  as.owin(X, ..., fatal=fatal)
}
