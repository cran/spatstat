#
#    as.im.R
#
#    conversion to class "im"
#
#    $Revision: 1.1 $   $Date: 2004/01/06 10:17:04 $
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
    return(im(v, w$xcol, w$yrow))
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
    if(!funnywindow) {
      values <- f(xx, yy, ...)
      v <- matrix(values, nrow=nrow(m), ncol=ncol(m))
    } else {
      xx <- xx[m]
      yy <- yy[m]
      values <- f(xx, yy, ...)
      v <- matrix(, nrow=nrow(m), ncol=ncol(m))
      v[m] <- values
      v[!m] <- NA
    }
    return(im(v, w$xcol, w$yrow))
  }

  stop("Can't convert x to a pixel image")
}
