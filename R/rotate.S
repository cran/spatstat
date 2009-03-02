#
#	rotate.S
#
#	$Revision: 1.12 $	$Date: 2009/02/27 18:47:37 $
#

rotxy <- function(X, angle=pi/2) {
  co <- cos(angle)
  si <- sin(angle)
  list(x = co * X$x - si * X$y,
       y = si * X$x + co * X$y)
}

rotxypolygon <- function(p, angle=pi/2) {
  p[c("x","y")] <- rotxy(p, angle=angle)
  # area and hole status are invariant under rotation
  return(p)
}

"rotate.owin" <- function(X, angle=pi/2, ...) {
  verifyclass(X, "owin")
  switch(X$type,
         rectangle={
           # convert rectangle to polygon
           P <- owin(X$xrange, X$yrange, poly=
                     list(x=X$xrange[c(1,2,2,1)],
                          y=X$yrange[c(1,1,2,2)]),
                     unitname=unitname(X))
           # call polygonal case
           return(rotate.owin(P, angle))
         },
         polygonal={
           # First rotate the polygonal boundaries
           bdry <- lapply(X$bdry, rotxypolygon, angle=angle)
           # wrap up
           W <- owin(poly=bdry, unitname=unitname(X))
           W <- rescue.rectangle(W)
           return(W)
         },
         mask={
           stop(paste("Sorry,", sQuote("rotate.owin"),
                      "is not yet implemented for masks"))
         },
         stop("Unrecognised window type")
         )
}

"rotate.ppp" <- function(X, angle=pi/2, ...) {
  verifyclass(X, "ppp")
  r <- rotxy(X, angle)
  w <- rotate.owin(X$window, angle)
  return(ppp(r$x, r$y, window=w, marks=marks(X, dfok=TRUE), check=FALSE))
}


"rotate" <- function(X, ...) {
  UseMethod("rotate")
}

  
