#
#  fryplot.R
#
#  $Revision: 1.4 $ $Date: 2010/11/26 09:07:11 $
#

fryplot <- function(X, ..., width=NULL) {
  Xname <- deparse(substitute(X))
  X <- as.ppp(X)
  b <- as.rectangle(X)
  halfspan <- with(b, c(diff(xrange), diff(yrange)))/2
  if(!is.null(width)) {
    halfwidth <- ensure2vector(width)/2
    halfspan <- pmin(halfspan, halfwidth)
  }
  bb <- owin(c(-1,1) * halfspan[1], c(-1,1) * halfspan[2])
  do.call("plot.owin",
          resolve.defaults(list(bb),
                           list(...),
                           list(invert=TRUE),
                           list(main=paste("Fry plot of", Xname))))
  n <- X$n
  xx <- X$x
  yy <- X$y
  for(i in 1:n) {
    dxi <- xx[-i] - xx[i]
    dyi <- yy[-i] - yy[i]
    oki <- (abs(dxi) < halfspan[1]) & (abs(dyi) < halfspan[2])
    if(any(oki)) 
      do.call.matched("points.default",
                      append(list(x=dxi[oki], y=dyi[oki]),
                             list(...)),
                      extrargs=c("pch", "col", "bg", "cex", "lwd"))
  }
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
