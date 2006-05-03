#
#    as.im.R
#
#    conversion to class "im"
#
#    $Revision: 1.6 $   $Date: 2006/04/11 10:20:21 $
#
#    as.im()
#
as.im <- function(X, W, ...) {

  x <- X
  
  if(verifyclass(x, "im", fatal=FALSE))
    return(x)

  if(verifyclass(x, "owin", fatal=FALSE)) {
    w <- as.mask(x)
    m <- w$m
    v <- m * 1
    v[!m] <- NA
    out <- list(v = v, 
                dim    = w$dim,
                xrange = w$xrange,
                yrange = w$yrange,
                xstep  = w$xstep,
                ystep  = w$ystep,
                xcol   = w$xcol,
                yrow   = w$yrow,
                lev    = NULL,
                type    = "integer")
    class(out) <- "im"
    return(out)
  }

  if(is.numeric(x) && length(x) == 1) {
    xvalue <- x
    x <- function(xx, yy, ...) { rep(xvalue, length(xx)) }
  }
  
  if(is.function(x)) {
    f <- x 
    w <- as.owin(W)
    w <- as.mask(w)
    m <- w$m
    funnywindow <- !all(m)

    xx <- raster.x(w)
    yy <- raster.y(w)
    lev <- NULL
    
    if(!funnywindow) {
      values <- f(xx, yy, ...)
      if(is.factor(values)) {
        lev <- levels(values)
        values <- as.integer(values)
      }
      v <- matrix(values, nrow=nrow(m), ncol=ncol(m))
    } else {
      xx <- xx[m]
      yy <- yy[m]
      values <- f(xx, yy, ...)
      if(is.factor(values)) {
        lev <- levels(values)
        values <- as.integer(values)
      }
      v <- matrix(, nrow=nrow(m), ncol=ncol(m))
      v[m] <- values
      v[!m] <- NA
    }
    return(im(v, w$xcol, w$yrow, lev))
  }

  if(is.list(x) && checkfields(x, c("x","y","z"))) {
    stopifnot(is.matrix(x$z))
    z <- x$z
    y <- x$y
    x <- x$x
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
    return(im(t(z), x, y))
  }
    
  stop("Can't convert x to a pixel image")
}
