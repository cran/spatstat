#
#  fryplot.R
#
#  $Revision: 1.2 $ $Date: 2008/12/09 22:13:31 $
#

fryplot <- function(X, ..., width=NULL) {
  Xname <- deparse(substitute(X))
  X <- as.ppp(X)
  b <- as.rectangle(X)
  bb <- owin(c(-1,1) * diff(b$xrange), c(-1,1) * diff(b$yrange))
  if(!is.null(width)){
    limits <- c(-1,1)*width/2
    xylim <- list(xlim=limits, ylim=limits)
  } else xylim <- NULL
  do.call("plot.owin",
          resolve.defaults(list(bb, type="n"),
                           list(...),
                           xylim,
                           list(main=paste("Fry plot of", Xname))))
  n <- X$n
  xx <- X$x
  yy <- X$y
  for(i in 1:n)
    do.call.matched("points.default",
                    append(list(x=xx[-i] - xx[i], y=yy[-i] - yy[i]),
                           list(...)),
                    extrargs=c("pch", "col", "bg", "cex", "lwd"))
  return(invisible(NULL))
}

frypoints <- function(X) {
  X <- as.ppp(X)
  b <- as.rectangle(X)
  bb <- owin(c(-1,1) * diff(b$xrange), c(-1,1) * diff(b$yrange))
  n <- X$n
  xx <- X$x
  yy <- X$y
  dx <- outer(xx, xx, "-")
  dy <- outer(yy, yy, "-")
  nondiag <- matrix(TRUE, n, n)
  diag(nondiag) <- FALSE
  DX <- as.vector(dx[nondiag])
  DY <- as.vector(dy[nondiag])
  Fry <- ppp(DX, DY, window=bb, check=FALSE)
  return(Fry)
}
