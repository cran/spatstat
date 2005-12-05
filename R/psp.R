#
#  psp.R
#
#  $Revision: 1.3 $ $Date: 2005/12/02 06:51:09 $
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

as.psp <- function(x, ...) {
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

#################################################
#  plot and print methods
#################################################

plot.psp <- function(x, ...) {
  verifyclass(x, "psp")
  do.call.matched("plot.owin", 
                  resolve.defaults(list(x=x$window),
                                   list(...),
                                   list(main=deparse(substitute(x)))))
  if(!is.null(x$marks))
    warning("marks are currently ignored in plot.psp")
  do.call.matched("segments", append(as.list(x$ends), list(...)))
}

print.psp <- function(x, ...) {
  verifyclass(x, "psp")
  cat(paste("planar line segment pattern:",
            x$n, "line segments\n"))
  print(x$window)
  if(!is.null(x$marks))
    cat(paste("Marks vector of type", sQuote(typeof(x$marks)), "\n"))
  return(invisible(NULL))
}

####################################################
#    summary information
####################################################

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
  function(x, subset, window, drop, ...) {

        verifyclass(x, "psp")

        trim <- !missing(window)
        thin <- !missing(subset)
        if(!thin && !trim)
          stop("Please specify a subset (to thin the pattern) or a window (to trim it)")

        # thin first, according to 'subset'
        if(thin)
          x <- as.psp(x$ends[subset,],
                       window=x$window,
                       marks=if(is.null(x$marks)) NULL else x$marks[subset])

        # now trim to window 
        if(trim) 
          x <- clip.psp(x, window=window)
        
        return(x)
}
  

########################################################
# clipping operation (for subset)
########################################################

clip.psp <- function(x, window, check=TRUE) {
  verifyclass(x, "psp")
  verifyclass(window, "owin")
  if(check && !is.subset.owin(window, x$window))
    warning("The clipping window is not a subset of the window containing the line segment pattern x")
  if(window$type != "rectangle")
    stop("sorry, clipping is only implemented for rectangular windows")
  ends <- x$ends
  marks <- x$marks
  # accept segments which are entirely inside the window
  # (by convexity)
  in0 <- inside.owin(ends$x0, ends$y0, window)
  in1 <- inside.owin(ends$x1, ends$y1, window)
  ok <- in0 & in1
  ends.inside <- ends[ok,]
  marks.inside <- if(is.null(marks)) NULL else marks[ok]
  # consider the rest
  ends <- ends[!ok,]
  in0 <- in0[!ok] 
  in1 <- in1[!ok]
  if(!is.null(marks)) marks <- marks[!ok]
  # first clip segments to the range x \in [xmin, xmax]
  # use parametric coordinates
  small <- function(x) { abs(x) <= .Machine$double.eps }
  tvalue <- function(z0, z1, zt) {
    y1 <- z1 - z0
    yt <- zt - z0
    tval <- ifelse(small(y1), 0.5, yt/y1)
    betwee <- (yt * (zt - z1)) <= 0
    return(ifelse(betwee, tval, NA))
  }
  between <- function(x, r) { ((x-r[1]) * (x-r[2])) <= 0 }
  tx <- cbind(ifelse(between(ends$x0, window$xrange), 0, NA),
              ifelse(between(ends$x1, window$xrange), 1, NA),
              tvalue(ends$x0, ends$x1, window$xrange[1]),
              tvalue(ends$x0, ends$x1, window$xrange[2]))
  # discard segments which do not lie in the x range 
  nx <- apply(!is.na(tx), 1, sum)
  ok <- (nx >= 2)
  ends <- ends[ok,]
  tx   <- tx[ok,]
  in0  <- in0[ok]
  in1  <- in1[ok]
  if(!is.null(marks)) marks <- marks[ok]
  # Clip the segments to the x range
  tmin <- apply(tx, 1, min, na.rm=TRUE)
  tmax <- apply(tx, 1, max, na.rm=TRUE)
  dx <- ends$x1 - ends$x0
  dy <- ends$y1 - ends$y0
  ends.xclipped <- data.frame(x0=ends$x0 + tmin * dx,
                             y0=ends$y0 + tmin * dy,
                             x1=ends$x0 + tmax * dx,
                             y1=ends$y0 + tmax * dy)
  # Now clip the segments to the range y \in [ymin, ymax]
  ends <- ends.xclipped
  in0 <- inside.owin(ends$x0, ends$y0, window)
  in1 <- inside.owin(ends$x1, ends$y1, window)
  ty <- cbind(ifelse(in0, 0, NA),
              ifelse(in1, 1, NA),
              tvalue(ends$y0, ends$y1, window$yrange[1]),
              tvalue(ends$y0, ends$y1, window$yrange[2]))
  # discard segments which do not lie in the y range 
  ny <- apply(!is.na(ty), 1, sum)
  ok <- (ny >= 2)
  ends <- ends[ok,]
  ty   <- ty[ok,]
  in0  <- in0[ok]
  in1  <- in1[ok]
  if(!is.null(marks)) marks <- marks[ok]
  # Clip the segments to the y range
  tmin <- apply(ty, 1, min, na.rm=TRUE)
  tmax <- apply(ty, 1, max, na.rm=TRUE)
  dx <- ends$x1 - ends$x0
  dy <- ends$y1 - ends$y0
  ends.clipped <- data.frame(x0=ends$x0 + tmin * dx,
                             y0=ends$y0 + tmin * dy,
                             x1=ends$x0 + tmax * dx,
                             y1=ends$y0 + tmax * dy)
  marks.clipped <- marks
  # OK - segments clipped
  # Put them together with the unclipped ones
  ends.all <- rbind(ends.inside, ends.clipped)
  marks.all <- c(marks.inside, marks.clipped)
  as.psp(ends.all, window=window, marks=marks)
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
  psp(ends0$x, ends0$y, ends1$x, ends1$y, window=W, marks=X$marks)
}

shift.psp <- function(X, ...) {
  verifyclass(X, "psp")
  W <- shift.owin(X$window, ...)
  E <- X$ends
  ends0 <- shiftxy(list(x=E$x0,y=E$y0), ...)
  ends1 <- shiftxy(list(x=E$x1,y=E$y1), ...)
  psp(ends0$x, ends0$y, ends1$x, ends1$y, window=W, marks=X$marks)
}

rotate.psp <- function(X, ...) {
  verifyclass(X, "psp")
  W <- rotate.owin(X$window, ...)
  E <- X$ends
  ends0 <- rotxy(list(x=E$x0,y=E$y0), ...)
  ends1 <- rotxy(list(x=E$x1,y=E$y1), ...)
  psp(ends0$x, ends0$y, ends1$x, ends1$y, window=W, marks=X$marks)
}


