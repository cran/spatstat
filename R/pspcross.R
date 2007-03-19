#
#    pspcross.R
#
#    Intersections of line segments
#    
#    $Revision: 1.3 $   $Date: 2007/03/16 13:32:41 $
#
#
crossing.psp <- function(A,B) {
  verifyclass(A, "psp")
  verifyclass(B, "psp")
  eps <- .Machine$double.eps

  na <- A$n
  eA <- A$ends
  x0a <- eA$x0
  y0a <- eA$y0
  dxa <- eA$x1 - eA$x0
  dya <- eA$y1 - eA$y0

  nb <- B$n
  eB <- B$ends
  x0b <- eB$x0
  y0b <- eB$y0
  dxb <- eB$x1 - eB$x0
  dyb <- eB$y1 - eB$y0
  
  out <- .C("xysegint",
            na=as.integer(na),
            x0a=as.double(x0a),
            y0a=as.double(y0a),
            dxa=as.double(dxa),
            dya=as.double(dya), 
            nb=as.integer(nb),
            x0b=as.double(x0b),
            y0b=as.double(y0b),
            dxb=as.double(dxb),
            dyb=as.double(dyb), 
            eps=as.double(eps),
            xx=as.double(numeric(na * nb)),
            yy=as.double(numeric(na * nb)),
            ta=as.double(numeric(na * nb)),
            tb=as.double(numeric(na * nb)),
            ok=as.integer(integer(na * nb)),
     PACKAGE="spatstat")

  ok <- (matrix(out$ok, na, nb) != 0)
  xx <- matrix(out$xx, na, nb)
  yy <- matrix(out$yy, na, nb)
  xx <- as.vector(xx[ok])
  yy <- as.vector(yy[ok])
  result <- ppp(xx, yy, window=intersect.owin(A$window, B$window), check=FALSE)
  return(result)
}

test.crossing.psp <- function(A,B) {
  verifyclass(A, "psp")
  verifyclass(B, "psp")
  eps <- .Machine$double.eps

  na <- A$n
  eA <- A$ends
  x0a <- eA$x0
  y0a <- eA$y0
  dxa <- eA$x1 - eA$x0
  dya <- eA$y1 - eA$y0

  nb <- B$n
  eB <- B$ends
  x0b <- eB$x0
  y0b <- eB$y0
  dxb <- eB$x1 - eB$x0
  dyb <- eB$y1 - eB$y0
  
  out <- .C("xysi",
            na=as.integer(na),
            x0a=as.double(x0a),
            y0a=as.double(y0a),
            dxa=as.double(dxa),
            dya=as.double(dya), 
            nb=as.integer(nb),
            x0b=as.double(x0b),
            y0b=as.double(y0b),
            dxb=as.double(dxb),
            dyb=as.double(dyb), 
            eps=as.double(eps),
            ok=as.integer(integer(1)),
     PACKAGE="spatstat")$ok
  return(out != 0)
}


selfcrossing.psp <- function(A) {
  verifyclass(A, "psp")
  eps <- .Machine$double.eps

  n <- A$n
  eA <- A$ends
  x0 <- eA$x0
  y0 <- eA$y0
  dx <- eA$x1 - eA$x0
  dy <- eA$y1 - eA$y0

  out <- .C("xysegXint",
            n=as.integer(n),
            x0=as.double(x0),
            y0=as.double(y0),
            dx=as.double(dx),
            dy=as.double(dy), 
            eps=as.double(eps),
            xx=as.double(numeric(n^2)),
            yy=as.double(numeric(n^2)),
            ti=as.double(numeric(n^2)),
            tj=as.double(numeric(n^2)),
            ok=as.integer(integer(n^2)),
     PACKAGE="spatstat")

  ok <- (matrix(out$ok, n, n) != 0)
  xx <- matrix(out$xx, n, n)
  yy <- matrix(out$yy, n, n)
  xx <- as.vector(xx[ok])
  yy <- as.vector(yy[ok])
  result <- ppp(xx, yy, window=A$window, check=FALSE)
  return(result)
}


