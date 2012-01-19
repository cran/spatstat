#
#   nncross.R
#
#
#    $Revision: 1.11 $  $Date: 2012/01/17 09:18:32 $
#


nncross <- function(X, Y, iX=NULL, iY=NULL) {
  X <- as.ppp(X, W=bounding.box.xy)
  stopifnot(is.ppp(Y) || is.psp(Y))

  # deal with null cases
  nX <- X$n
  nY <- Y$n
  if(nX == 0)
    return(data.frame(dist=numeric(0), which=integer(0)))
  if(nY == 0)
    return(data.frame(dist=rep(Inf, nX), which=rep(NA, nX)))

  # Y is a line segment pattern 
  if(is.psp(Y))
    return(ppllengine(X,Y,"distance"))

  # Y is a point pattern
  if(is.null(iX) != is.null(iY))
    stop("If one of iX, iY is given, then both must be given")
  exclude <- (!is.null(iX) || !is.null(iY))
  if(exclude) {
    stopifnot(is.integer(iX) && is.integer(iY))
    if(length(iX) != nX)
      stop("length of iX does not match the number of points in X")
    if(length(iY) != nY)
      stop("length of iY does not match the number of points in Y")
  }
    
  # sort in increasing order of y coordinate
  oX <- order(X$y)
  X <- X[oX]
  oY <- order(Y$y)
  Y <- Y[oY]
  if(exclude) {
    iX <- iX[oX]
    iY <- iY[oY]
  }

  # call C code
  nndv <- numeric(nX)
  nnwh <- integer(nX)

  DUP <- spatstat.options("dupC")
  huge <- 1.1 * diameter(bounding.box(as.rectangle(X), as.rectangle(Y)))
  
  if(!exclude) 
    z <- .C("nnXwhich",
            n1=as.integer(nX),
            x1=as.double(X$x),
            y1=as.double(X$y),
            n2=as.integer(nY),
            x2=as.double(Y$x),
            y2=as.double(Y$y),
            nnd=as.double(nndv),
            nnwhich=as.integer(nnwh),
            huge=as.double(huge),
            DUP=DUP,
            PACKAGE="spatstat")
  else
    z <- .C("nnXexclude",
            n1=as.integer(nX),
            x1=as.double(X$x),
            y1=as.double(X$y),
            id1=as.integer(iX),
            n2=as.integer(nY),
            x2=as.double(Y$x),
            y2=as.double(Y$y),
            id2=as.integer(iY),
            nnd=as.double(nndv),
            nnwhich=as.integer(nnwh),
            huge=as.double(huge),
            DUP=DUP,
            PACKAGE="spatstat")
    
  # reinterpret in original ordering
  nndv[oX] <- z$nnd
  nnwcode <- z$nnwhich + 1
  if(any(uhoh <- (nnwcode < 1))) {
    warning("NA's produced in nncross()$which")
    nnwcode[uhoh] <- NA
  }
  nnwh[oX] <- oY[nnwcode]
  return(data.frame(dist=nndv, which=nnwh))
}

