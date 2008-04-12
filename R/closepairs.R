#
# closepairs.R
#
#   $Revision: 1.7 $   $Date: 2008/04/02 13:38:54 $
#
#  simply extract the r-close pairs from a dataset
# 
#  Less memory-hungry for large patterns
#
closepairs <- function(X, rmax) {
  verifyclass(X, "ppp")
  stopifnot(is.numeric(rmax) && length(rmax) == 1 && rmax >= 0)
  oo <- order(X$x)
  Xsort <- X[oo]
  DUP <- spatstat.options("dupC")
  npairs <- .C("paircount",
            nxy=as.integer(X$n),
            x=as.double(Xsort$x),
            y=as.double(Xsort$y),
            rmaxi=as.double(rmax),
            count=as.integer(integer(1)),
            DUP=DUP,
            PACKAGE="spatstat")$count
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

  z <-
    .C("closepairs",
       nxy=as.integer(X$n),
       x=as.double(Xsort$x),
       y=as.double(Xsort$y),
       r=as.double(rmax),
       noutmax=as.integer(npairs), 
       nout=as.integer(integer(1)),
       iout=as.integer(integer(npairs)),
       jout=as.integer(integer(npairs)), 
       xiout=as.double(numeric(npairs)),
       yiout=as.double(numeric(npairs)),
       xjout=as.double(numeric(npairs)),
       yjout=as.double(numeric(npairs)),
       dxout=as.double(numeric(npairs)),
       dyout=as.double(numeric(npairs)),
       dout=as.double(numeric(npairs)),
       status=as.integer(integer(1)),
       DUP=DUP,
       PACKAGE="spatstat")
  if(z$status != 0)
    stop(paste("Internal error: C routine complains that insufficient space was allocated:", npairs))
  if(z$nout != npairs)
    warning(paste("Internal error: npairs miscounted:", npairs, "!= ", z$nout))
  # trim vectors to the length indicated
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
  ooX <- order(X$x)
  Xsort <- X[ooX]
  ooY <- order(Y$x)
  Ysort <- Y[ooY]
  DUP <- spatstat.options("dupC")
  npairs <- .C("crosscount",
               nn1=as.integer(X$n),
               x1=as.double(Xsort$x),
               y1=as.double(Xsort$y),
               nn2=as.integer(Ysort$n),
               x2=as.double(Ysort$x),
               y2=as.double(Ysort$y),
               rmaxi=as.double(rmax),
               count=as.integer(integer(1)),
               DUP=DUP,
               PACKAGE="spatstat")$count
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

  z <-
    .C("crosspairs",
       nn1=as.integer(X$n),
       x1=as.double(Xsort$x),
       y1=as.double(Xsort$y),
       nn2=as.integer(Y$n),
       x2=as.double(Ysort$x),
       y2=as.double(Ysort$y),
       r=as.double(rmax),
       noutmax=as.integer(npairs), 
       nout=as.integer(integer(1)),
       iout=as.integer(integer(npairs)),
       jout=as.integer(integer(npairs)), 
       xiout=as.double(numeric(npairs)),
       yiout=as.double(numeric(npairs)),
       xjout=as.double(numeric(npairs)),
       yjout=as.double(numeric(npairs)),
       dxout=as.double(numeric(npairs)),
       dyout=as.double(numeric(npairs)),
       dout=as.double(numeric(npairs)),
       status=as.integer(integer(1)),
       DUP=DUP,
       PACKAGE="spatstat")
  if(z$status != 0)
    stop(paste("Internal error: C routine complains that insufficient space was allocated:", npairs))
  if(z$nout != npairs)
    warning(paste("Internal error: npairs miscounted:", npairs, "!= ", z$nout))
  # trim vectors to the length indicated
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
