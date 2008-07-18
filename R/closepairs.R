#
# closepairs.R
#
#   $Revision: 1.11 $   $Date: 2008/07/02 00:51:35 $
#
#  simply extract the r-close pairs from a dataset
# 
#  Less memory-hungry for large patterns
#
closepairs <- function(X, rmax) {
  verifyclass(X, "ppp")
  stopifnot(is.numeric(rmax) && length(rmax) == 1 && rmax >= 0)
  npoints <- X$n
  if(npoints == 0)
    return(list(i=integer(0),
                j=integer(0),
                xi=numeric(0),
                yi=numeric(0),
                xj=numeric(0),
                yj=numeric(0),
                dx=numeric(0),
                dy=numeric(0),
                d=numeric(0)))
  # sort points by increasing x coordinate
  oo <- order(X$x)
  Xsort <- X[oo]
  # First make an OVERESTIMATE of the number of pairs
  nsize <- ceiling(4 * pi * (npoints^2) * (rmax^2)/area.owin(X$window))
  nsize <- max(1024, nsize)
  # Now hopefully extract pairs
  DUP <- spatstat.options("dupC")
  z <-
    .C("closepairs",
       nxy=as.integer(npoints),
       x=as.double(Xsort$x),
       y=as.double(Xsort$y),
       r=as.double(rmax),
       noutmax=as.integer(nsize), 
       nout=as.integer(integer(1)),
       iout=as.integer(integer(nsize)),
       jout=as.integer(integer(nsize)), 
       xiout=as.double(numeric(nsize)),
       yiout=as.double(numeric(nsize)),
       xjout=as.double(numeric(nsize)),
       yjout=as.double(numeric(nsize)),
       dxout=as.double(numeric(nsize)),
       dyout=as.double(numeric(nsize)),
       dout=as.double(numeric(nsize)),
       status=as.integer(integer(1)),
       DUP=DUP,
       PACKAGE="spatstat")

  if(z$status != 0) {
    # Guess was insufficient
    # Obtain an OVERCOUNT of the number of pairs
    # (to work around gcc bug #323)
    rmaxplus <- 1.25 * rmax
    nsize <- .C("paircount",
            nxy=as.integer(npoints),
            x=as.double(Xsort$x),
            y=as.double(Xsort$y),
            rmaxi=as.double(rmaxplus),
            count=as.integer(integer(1)),
            DUP=DUP,
            PACKAGE="spatstat")$count
    if(nsize <= 0)
      return(list(i=integer(0),
                  j=integer(0),
                  xi=numeric(0),
                  yi=numeric(0),
                  xj=numeric(0),
                  yj=numeric(0),
                  dx=numeric(0),
                  dy=numeric(0),
                  d=numeric(0)))
    # add a bit more for safety
    nsize <- ceiling(1.1 * nsize) + 2 * npoints
    # now extract points
    z <-
      .C("closepairs",
         nxy=as.integer(npoints),
         x=as.double(Xsort$x),
         y=as.double(Xsort$y),
         r=as.double(rmax),
         noutmax=as.integer(nsize), 
         nout=as.integer(integer(1)),
         iout=as.integer(integer(nsize)),
         jout=as.integer(integer(nsize)), 
         xiout=as.double(numeric(nsize)),
         yiout=as.double(numeric(nsize)),
         xjout=as.double(numeric(nsize)),
         yjout=as.double(numeric(nsize)),
         dxout=as.double(numeric(nsize)),
         dyout=as.double(numeric(nsize)),
         dout=as.double(numeric(nsize)),
         status=as.integer(integer(1)),
         DUP=DUP,
         PACKAGE="spatstat")
    if(z$status != 0)
      stop(paste("Internal error: C routine complains that insufficient space was allocated:", nsize))
  }
  
  # trim vectors to the length indicated
  npairs <- z$nout
  if(npairs <= 0)
    return(list(i=integer(0),
                j=integer(0),
                xi=numeric(0),
                yi=numeric(0),
                xj=numeric(0),
                yj=numeric(0),
                dx=numeric(0),
                dy=numeric(0),
                d=numeric(0)))
  actual <- seq(npairs)
  i  <- z$iout[actual] + 1
  j  <- z$jout[actual] + 1
  xi <- z$xiout[actual]
  yi <- z$yiout[actual]
  xj <- z$xjout[actual]
  yj <- z$yjout[actual]
  dx <- z$dxout[actual]
  dy <- z$dyout[actual]
  d <-  z$dout[actual]
  # convert i,j indices to original sequence
  i <- oo[i]
  j <- oo[j]
  # done
  return(list(i=i,
              j=j,
              xi=xi, 
              yi=yi,
              xj=xj,
              yj=yj,
              dx=dx,
              dy=dy,
              d=d))
}


#######################

crosspairs <- function(X, Y, rmax) {
  verifyclass(X, "ppp")
  verifyclass(Y, "ppp")
  stopifnot(is.numeric(rmax) && length(rmax) == 1 && rmax >= 0)
  # order patterns by increasing x coordinate
  ooX <- order(X$x)
  Xsort <- X[ooX]
  ooY <- order(Y$x)
  Ysort <- Y[ooY]
  # obtain upper estimate of number of pairs
  # (to work around gcc bug 323)
  DUP <- spatstat.options("dupC")
  rmaxplus <- 1.25 * rmax
  nsize <- .C("crosscount",
               nn1=as.integer(X$n),
               x1=as.double(Xsort$x),
               y1=as.double(Xsort$y),
               nn2=as.integer(Ysort$n),
               x2=as.double(Ysort$x),
               y2=as.double(Ysort$y),
               rmaxi=as.double(rmaxplus),
               count=as.integer(integer(1)),
               DUP=DUP,
               PACKAGE="spatstat")$count
  if(nsize <= 0)
    return(list(i=integer(0),
                j=integer(0),
                xi=numeric(0),
                yi=numeric(0),
                xj=numeric(0),
                yj=numeric(0),
                dx=numeric(0),
                dy=numeric(0),
                d=numeric(0)))

  # allow slightly more space to work around gcc bug #323
  nsize <- ceiling(1.1 * nsize) + X$n + Y$n

  # now extract pairs
  z <-
    .C("crosspairs",
       nn1=as.integer(X$n),
       x1=as.double(Xsort$x),
       y1=as.double(Xsort$y),
       nn2=as.integer(Y$n),
       x2=as.double(Ysort$x),
       y2=as.double(Ysort$y),
       r=as.double(rmax),
       noutmax=as.integer(nsize), 
       nout=as.integer(integer(1)),
       iout=as.integer(integer(nsize)),
       jout=as.integer(integer(nsize)), 
       xiout=as.double(numeric(nsize)),
       yiout=as.double(numeric(nsize)),
       xjout=as.double(numeric(nsize)),
       yjout=as.double(numeric(nsize)),
       dxout=as.double(numeric(nsize)),
       dyout=as.double(numeric(nsize)),
       dout=as.double(numeric(nsize)),
       status=as.integer(integer(1)),
       DUP=DUP,
       PACKAGE="spatstat")
  if(z$status != 0)
    stop(paste("Internal error: C routine complains that insufficient space was allocated:", nsize))
  # trim vectors to the length indicated
  npairs <- z$nout
  if(npairs <= 0)
    return(list(i=integer(0),
                j=integer(0),
                xi=numeric(0),
                yi=numeric(0),
                xj=numeric(0),
                yj=numeric(0),
                dx=numeric(0),
                dy=numeric(0),
                d=numeric(0)))
  actual <- seq(npairs)
  i  <- z$iout[actual] + 1
  j  <- z$jout[actual] + 1
  xi <- z$xiout[actual]
  yi <- z$yiout[actual]
  xj <- z$xjout[actual]
  yj <- z$yjout[actual]
  dx <- z$dxout[actual]
  dy <- z$dyout[actual]
  d <-  z$dout[actual]
  # convert i,j indices to original sequences
  i <- ooX[i]
  j <- ooY[j]
  # done
  return(list(i=i,
              j=j,
              xi=xi, 
              yi=yi,
              xj=xj,
              yj=yj,
              dx=dx,
              dy=dy,
              d=d))
}
