#
#  psp.R
#
#  $Revision: 1.23 $ $Date: 2006/10/13 09:38:04 $
#
# Class "psp" of planar line segment patterns
#
#
#################################################
# creator
#################################################
psp <- function(x0, y0, x1, y1, window, marks=NULL) {
  stopifnot(is.numeric(x0))
  stopifnot(is.numeric(y0))
  stopifnot(is.numeric(x1))
  stopifnot(is.numeric(y1))
  stopifnot(is.vector(x0))
  stopifnot(is.vector(y0))
  stopifnot(is.vector(x1))
  stopifnot(is.vector(y1))
  stopifnot(length(x0) == length(y0))
  stopifnot(length(x1) == length(y1))
  stopifnot(length(x0) == length(x1))
  ends <- data.frame(x0=x0,y0=y0,x1=x1,y1=y1)
  window <- as.owin(window)
  if(!is.null(marks)) {
    if(is.data.frame(marks))
      stop("Sorry, data frames of marks are not implemented for psp objects")
    stopifnot(is.vector(marks))
    stopifnot(length(marks) == length(x0))
  }
  out <- list(ends=ends,
              window=window,
              n = nrow(ends),
              marks = marks)
  class(out) <- c("psp", class(out))
  return(out)
}

######################################################
#  conversion
######################################################

as.psp <- function(x, ..., from=NULL, to=NULL) {
  # special case: two point patterns
  if(is.null(from) != is.null(to))
    stop(paste("If one of", sQuote("from"), "and", sQuote("to"),
               "is specified, then both must be specified"))
  if(!is.null(from) && !is.null(to)) {
    verifyclass(from, "ppp")
    verifyclass(to, "ppp")
    if(from$n != to$n)
      stop(paste("point patterns", sQuote("from"), "and", sQuote("to"),
                 "have different numbers of points"))
    uni <- union.owin(from$window, to$window)
    Y <- do.call("psp",
                 resolve.defaults(list(from$x, from$y, to$x, to$y),
                                  list(...),
                                  list(window=uni)))
    return(Y)
  }
  UseMethod("as.psp")
}

as.psp.psp <- function(x, ..., fatal=TRUE) {
  if(!verifyclass(x, "psp", fatal=fatal))
    return(NULL)
  ends <- x$ends
  psp(ends$x0, ends$y0, ends$x1, ends$y1, window=x$window, marks=x$marks)
}

as.psp.data.frame <- function(x, ..., window=NULL, marks=NULL, fatal=TRUE) {
  if(checkfields(x, c("x0", "y0", "x1", "y1")))
    return(psp(x$x0, x$y0, x$x1, x$y1, window=window, marks=marks))
  else if(checkfields(x, c("xmid", "ymid", "length", "angle"))) {
    dx <- cos(x$angle) * x$length/2
    dy <- sin(x$angle) * x$length/2
    return(psp(x$x - dx, x$y - dy, x$x + dx, x$y + dy,
                window=window, marks=marks))
  }
  else if(ncol(x) == 4)
    return(psp(x[,1], x[,2], x[,3], x[,4], window=window, marks=marks))
  else if(fatal)
    stop("Unable to interpret x as a line segment pattern")
  return(NULL)
}

as.psp.matrix <- function(x, ..., window=NULL, marks=NULL, fatal=TRUE) {
  if(ncol(x) == 4)
    return(psp(x[,1], x[,2], x[,3], x[,4], window=window, marks=marks))
  else if(fatal)
    stop("Unable to interpret x as a line segment pattern")
  return(NULL)
}

as.psp.default <- function(x, ..., window=NULL, marks=NULL, fatal=TRUE) {
  if(checkfields(x, c("x0", "y0", "x1", "y1")))
    return(psp(x$x0, x$y0, x$x1, x$y1, window=window, marks=marks))
  else if(checkfields(x, c("xmid", "ymid", "length", "angle"))) {
    dx <- cos(x$angle) * x$length/2
    dy <- sin(x$angle) * x$length/2
    return(psp(x$x - dx, x$y - dy, x$x + dx, x$y + dy,
                window=window, marks=marks))
  }
  else if(fatal)
    stop("Unable to interpret x as a line segment pattern")
  return(NULL)
}

as.psp.owin <- function(x, ..., fatal=TRUE) {
  verifyclass(x, "owin")
  # can't use as.rectangle here; still testing validity
  xframe <- owin(x$xrange, x$yrange)
  switch(x$type,
         rectangle = {
           xx <- x$xrange[c(1,2,2,1)]
           yy <- x$yrange[c(1,1,2,2)]
           nxt <- c(2,3,4,1)
           out <- psp(xx, yy, xx[nxt], yy[nxt], window=x)
           return(out)
         },
         polygonal = {
           x0 <- y0 <- x1 <- y1 <- numeric(0)
           bdry <- x$bdry
           for(i in seq(bdry)) {
             po <- bdry[[i]]
             ni <- length(po$x)
             nxt <- c(2:ni, 1)
             x0 <- c(x0, po$x)
             y0 <- c(y0, po$y)
             x1 <- c(x1, po$x[nxt])
             y1 <- c(y1, po$y[nxt])
           }
           out <- psp(x0, y0, x1, y1,  window=xframe)
           return(out)
         },
         mask = {
           if(fatal) stop("x is a mask")
           else warning("x is a mask - no line segments returned")
           return(psp(numeric(0), numeric(0), numeric(0), numeric(0),
                      window=xframe))
         })
  return(NULL)
}


append.psp <- function(A,B) {
  verifyclass(A, "psp")
  verifyclass(B, "psp")
  stopifnot(identical(A$window, B$window))
  result <- list(ends=rbind(A$ends, B$ends),
                 window=A$window,
                 n=A$n + B$n,
                 marks=c(marks(A), marks(B)))
  class(result) <- c("psp", class(result))
  return(result)
}

rebound.psp <- function(x, ...) {
  verifyclass(x, "psp")
  x$window <- rebound.owin(x$window, ...)
  return(x)
}


marks.psp <- function(x, ..., dfok=FALSE) {
  # data frames of marks are not implemented for psp
  return(x$marks)
}

markformat.psp <- function(x) {
  # data frames of marks are not implemented for psp
  return(if(!is.null(marks(x))) "vector" else "none")
}

#################################################
#  plot and print methods
#################################################

plot.psp <- function(x, ..., add=FALSE) {
  verifyclass(x, "psp")
  if(!add) {
    do.call.matched("plot.owin", 
                    resolve.defaults(list(x=x$window),
                                     list(...),
                                     list(main=deparse(substitute(x)))))
  }
  if(!is.null(marks(x, dfok=TRUE)))
    warning("marks are currently ignored in plot.psp")
  if(x$n > 0)
    do.call.matched("segments", append(as.list(x$ends), list(...)))
  return(invisible(NULL))
}

print.psp <- function(x, ...) {
  verifyclass(x, "psp")
  cat(paste("planar line segment pattern:",
            x$n, "line segments\n"))
  print(x$window)
  marx <- marks(x, dfok=FALSE)
  if(!is.null(marx))
    cat(paste("Marks vector of type", sQuote(typeof(marx)), "\n"))
  return(invisible(NULL))
}

####################################################
#    summary information
####################################################

endpoints.psp <- function(x, which="both") {
  verifyclass(x, "psp")
  ends <- x$ends
  n <- x$n
  switch(which,
         both={
           first <- second <- rep(TRUE, n)
         },
         first={
           first <- rep(TRUE, n)
           second <- rep(FALSE, n)
         },
         second={
           first <- rep(FALSE, n)
           second <- rep(TRUE, n)
         },
         left={
           first <- (ends$x0 < ends$x1)
           second <- !first
         },
         right={
           first <- (ends$x0 > ends$x1)
           second <- !first
         },
         lower={
           first <- (ends$y0 < ends$y1)
           second <- !first
         },
         upper={
           first <- (ends$y0 > ends$y1)
           second <- !first
         },
         stop(paste("Unrecognised option: which=", sQuote(which)))
         )
  ok <- rbind(first, second)
  xmat <- rbind(ends$x0, ends$x1)
  ymat <- rbind(ends$y0, ends$y1)
  idmat <- col(ok)
  xx <- as.vector(xmat[ok])
  yy <- as.vector(ymat[ok])
  id <- as.vector(idmat[ok])
  result <- ppp(xx, yy, window=x$window)
  attr(result, "id") <- id
  return(result)
}

midpoints.psp <- function(x) {
  verifyclass(x, "psp")
  xm <- eval(expression((x0+x1)/2), envir=x$ends)
  ym <- eval(expression((y0+y1)/2), envir=x$ends)
  ppp(x=xm, y=ym, window=x$window)
}

lengths.psp <- function(x) {
  verifyclass(x, "psp")
  eval(expression(sqrt((x1-x0)^2 + (y1-y0)^2)), envir=x$ends)
}

angles.psp <- function(x, directed=FALSE) {
  verifyclass(x, "psp")
  a <- eval(expression(atan2(y1-y0, x1-x0)), envir=x$ends)
  if(!directed) 
    a <- a %% pi
  return(a)
}

summary.psp <- function(object, ...) {
  verifyclass(object, "psp")
  len <- lengths.psp(object)
  out <- list(n = object$n,
              len = summary(len),
              totlen = sum(len),
              ang= summary(angles.psp(object)),
              w = summary.owin(object$window),
              marks=if(is.null(object$marks)) NULL else summary(object$marks))
  class(out) <- c("summary.psp", class(out))
  return(out)
}

print.summary.psp <- function(x, ...) {
  cat(paste(x$n, "line segments\n"))
  cat("Lengths:\n")
  print(x$len)
  cat(paste("Total length:", x$totlen, "\n"))
  cat(paste("Length per unit area:", x$totlen/x$w$area, "\n"))
  cat("Angles (radians):\n")
  print(x$ang)
  print(x$w)
  if(!is.null(x$marks)) {
    cat("Marks:\n")
    print(x$marks)
  }
  return(invisible(NULL))
}

  
########################################################
#  subsets
########################################################

"[.psp" <-
"subset.psp" <-
  function(x, i, j, drop, ...) {

    verifyclass(x, "psp")
    
    if(missing(i) && missing(j))
      return(x)
        
    if(!missing(i)) {
      style <- if(inherits(i, "owin")) "window" else "index"
      switch(style,
             window={
               x <- clip.psp(x, window=i)
             },
             index={
               subset <- i
               marktype <- markformat(x)
               x <- as.psp(x$ends[subset, ],
                           window=x$window,
                           marks=switch(marktype,
                             none=NULL,
                             vector=x$marks[subset],
                             dataframe=x$marks[subset,]))
             })
    }

    if(!missing(j))
      x <- x[j] # invokes code above
    
    return(x)
 }
  


#################################################
#  distance map
#################################################

distmap.psp <- function(X, ...) {
  verifyclass(X, "psp")
  W <- as.mask(X$window, ...)
  P <- cbind(as.vector(raster.x(W)), as.vector(raster.y(W)))
  U <- distppll(P, X$ends, mintype=2)
  xc <- W$xcol
  yr <- W$yrow
  Dist <- im(array(U$min.d, dim=W$dim), xc, yr)
  Indx <- im(array(U$min.which, dim=W$dim), xc, yr)
  Bdry <- im(bdist.pixels(W, coords=FALSE), xc, yr)
  attr(Dist, "index") <- Indx
  attr(Dist, "bdry")  <- Bdry
  return(Dist)
}

####################################################
# affine transformations
####################################################

affine.psp <- function(X,  ...) {
  verifyclass(X, "psp")
  W <- affine.owin(X$window, ...)
  E <- X$ends
  ends0 <- affinexy(list(x=E$x0,y=E$y0), ...)
  ends1 <- affinexy(list(x=E$x1,y=E$y1), ...)
  psp(ends0$x, ends0$y, ends1$x, ends1$y, window=W, marks=marks(X, dfok=TRUE))
}

shift.psp <- function(X, ...) {
  verifyclass(X, "psp")
  W <- shift.owin(X$window, ...)
  E <- X$ends
  ends0 <- shiftxy(list(x=E$x0,y=E$y0), ...)
  ends1 <- shiftxy(list(x=E$x1,y=E$y1), ...)
  psp(ends0$x, ends0$y, ends1$x, ends1$y, window=W, marks=marks(X, dfok=TRUE))
}

rotate.psp <- function(X, ...) {
  verifyclass(X, "psp")
  W <- rotate.owin(X$window, ...)
  E <- X$ends
  ends0 <- rotxy(list(x=E$x0,y=E$y0), ...)
  ends1 <- rotxy(list(x=E$x1,y=E$y1), ...)
  psp(ends0$x, ends0$y, ends1$x, ends1$y, window=W, marks=marks(X, dfok=TRUE))
}


