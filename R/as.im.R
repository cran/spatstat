#
#    as.im.R
#
#    conversion to class "im"
#
#    $Revision: 1.13 $   $Date: 2008/09/30 18:33:20 $
#
#    as.im()
#

as.im <- function(X, ...) {
  UseMethod("as.im")
}

as.im.im <- function(X,
                     W=as.mask(as.owin(X), dimyx=dimyx), ...,
                     dimyx=NULL, na.replace=NULL) {

  if(missing(W) && is.null(dimyx))
    return(na.handle.im(X, na.replace))

  # reshape pixel raster
  # invoke W = as.mask(X, dimyx)
  Y <- as.im(W, dimyx=dimyx)
  phase <- c((Y$xcol[1] - X$xcol[1])/X$xstep,
             (Y$yrow[1] - X$yrow[1])/X$ystep)
  Y$v <- matrixsample(X$v, Y$dim, phase=round(phase))

  return(na.handle.im(Y, na.replace))
}

as.im.owin <- function(X,
                     W=as.mask(as.owin(X), dimyx=dimyx), ...,
                     dimyx=NULL, na.replace=NULL) {

    # if W is missing, the default is now evaluated, as above.
    # if W is present, it may have to be converted
  if(!missing(W)) {
    stopifnot(is.owin(W))
    if(W$type != "mask")
      W <- as.mask(W, dimyx=dimyx)
  }
  m <- W$m
  v <- m * 1
  v[!m] <- NA
  out <- list(v = v, 
              dim    = W$dim,
              xrange = W$xrange,
              yrange = W$yrange,
              xstep  = W$xstep,
              ystep  = W$ystep,
              xcol   = W$xcol,
              yrow   = W$yrow,
              lev    = NULL,
              type    = "integer",
              units  = unitname(X))
  class(out) <- "im"
  return(na.handle.im(out, na.replace))
}


as.im.function <- function(X,
                          W=as.mask(as.owin(X), dimyx=dimyx), ...,
                          dimyx=NULL, na.replace=NULL) {
  f <- X
  W <- as.owin(W)
  W <- as.mask(W, dimyx=dimyx)
  m <- W$m
  funnywindow <- !all(m)

  xx <- as.vector(raster.x(W))
  yy <- as.vector(raster.y(W))
  lev <- NULL

  # evaluate function value at each pixel 
  if(!funnywindow) 
    values <- f(xx, yy, ...)
  else {
    # evaluate only inside window
    inside <- as.vector(m)
    val <- f(xx[inside], yy[inside], ...)
    # create space for full matrix
    msize <- length(m)
    values <-
      if(!is.factor(val))
        vector(mode=typeof(val), length=msize)
      else {
        lev <- levels(val)
        factor(rep(lev[1], msize), levels=lev)
      }
    # copy values, assigning NA outside window
    values[inside] <- val
    values[!inside] <- NA
  }
  
  if(is.factor(values)) 
    lev <- levels(values)
  
  out <- im(values, W$xcol, W$yrow, lev, unitname=unitname(W))
  return(na.handle.im(out, na.replace))
}

as.im.default <- function(X,
                          W=as.mask(as.owin(X), dimyx=dimyx), ...,
                          dimyx=NULL, na.replace=NULL) {

  if((is.vector(X) || is.factor(X)) && length(X) == 1) {
    # numerical value: interpret as constant function
    xvalue <- X
    X <- function(xx, yy, ...) { rep(xvalue, length(xx)) }
    return(as.im(X, W, ..., dimyx=dimyx, na.replace=na.replace))
  }
  
  if(is.list(X) && checkfields(X, c("x","y","z"))) {
    stopifnot(is.matrix(X$z))
    z <- X$z
    y <- X$y
    x <- X$x
    # Usual S convention as in contour.default() and image.default()
    # Rows of z correspond to x values.
    nr <- nrow(z)
    nc <- ncol(z)
    lx <- length(x)
    ly <- length(y)
    if(lx == nr + 1)
      x <- (x[-1] + x[-lx])/2
    else if(lx != nr)
      stop("length of x coordinate vector does not match number of rows of z")
    if(ly == nc + 1)
      y <- (y[-1] + y[-ly])/2
    else if(ly != nc)
      stop("length of y coordinate vector does not match number of columns of z")
    # convert to class "im"
    out <- im(t(z), x, y)
    # now apply W and dimyx if present
    if(missing(W) && !is.null(dimyx))
      out <- as.im(out, dimyx=dimyx)
    else if(!missing(W))
      out <- as.im(out, W=W, dimyx=dimyx)
    return(na.handle.im(out, na.replace))
  }
  stop("Can't convert X to a pixel image")
}

na.handle.im <- function(X, na.replace) {
if(is.null(na.replace))
  return(X)
if(length(na.replace) != 1)
  stop("na.replace should be a single value")
X$v[is.na(X$v)] <- na.replace
return(X)
}
