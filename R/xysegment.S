#
#      xysegment.S
#
#     $Revision: 1.2 $    $Date: 2001/11/19 05:08:06 $
#
# Low level utilities for analytic geometry for line segments
#
# author: Adrian Baddeley 2001
#         from an original by Rob Foxall 1997
#
# distpl(p, l) 
#       distance from a single point p  = (xp, yp)
#       to a single line segment l = (x1, y1, x2, y2)
#
# distppl(p, l) 
#       distances from each of a list of points p[i,]
#       to a single line segment l = (x1, y1, x2, y2)
#       [uses only vector parallel ops]
#
# distppll(p, l) 
#       distances from each of a list of points p[i,]
#       to each of a list of line segments l[i,] 
#       [uses large matrices and 'outer()']
#

distpl <- function(p, l) {
  xp <- p[1]
  yp <- p[2]
  dx <- l[3]-l[1]
  dy <- l[4]-l[2]
  leng <- sqrt(dx^2 + dy^2)
  # vector from 1st endpoint to p
  xpl <- xp - l[1]
  ypl <- yp - l[2]
  # distance from p to 1st & 2nd endpoints
  d1 <- sqrt(xpl^2 + ypl^2)
  d2 <- sqrt((xp-l[3])^2 + (yp-l[4])^2)
  dmin <- min(d1,d2)
  # rotation sine & cosine
  co <- dx/leng
  si <- dy/leng
  # back-rotated coords of p
  xpr <- co * xpl + si * ypl
  ypr <-  - si * xpl + co * ypl
  # test
  if(xpr >= 0 && xpr <= leng)
    dmin <- min(dmin, abs(ypr))
  return(dmin)
}

distppl <- function(p, l) {
  xp <- p[,1]
  yp <- p[,2]
  dx <- l[3]-l[1]
  dy <- l[4]-l[2]
  leng <- sqrt(dx^2 + dy^2)
  # vector from 1st endpoint to p
  xpl <- xp - l[1]
  ypl <- yp - l[2]
  # distance from p to 1st & 2nd endpoints
  d1 <- sqrt(xpl^2 + ypl^2)
  d2 <- sqrt((xp-l[3])^2 + (yp-l[4])^2)
  dmin <- pmin(d1,d2)
  # rotation sine & cosine
  co <- dx/leng
  si <- dy/leng
  # back-rotated coords of p
  xpr <- co * xpl + si * ypl
  ypr <-  - si * xpl + co * ypl
  # ypr is perpendicular distance to infinite line
  # Applies only when xp, yp in the middle
  middle <- (xpr >= 0 & xpr <= leng)
  if(any(middle))
    dmin[middle] <- pmin(dmin[middle], abs(ypr[middle]))
  
  return(dmin)
}

distppll <- function(p, l) {
  np <- nrow(p)
  nl <- nrow(l)
  xp <- p[,1]
  yp <- p[,2]
  dx <- l[,3]-l[,1]
  dy <- l[,4]-l[,2]
  # segment lengths
  leng <- sqrt(dx^2 + dy^2)
  # rotation sines & cosines
  co <- dx/leng
  si <- dy/leng
  co <- matrix(co, nrow=np, ncol=nl, byrow=T)
  si <- matrix(si, nrow=np, ncol=nl, byrow=T)
  # matrix of squared distances from p[i] to 1st endpoint of segment j
  xp.x1 <- outer(xp, l[,1], "-")
  yp.y1 <- outer(yp, l[,2], "-")
  d1 <- xp.x1^2 + yp.y1^2
  # ditto for 2nd endpoint
  xp.x2 <- outer(xp, l[,3], "-")
  yp.y2 <- outer(yp, l[,4], "-")
  d2 <- xp.x2^2 + yp.y2^2
  # for each (i,j) rotate p[i] around 1st endpoint of segment j
  # so that line segment coincides with x axis
  xpr <- xp.x1 * co + yp.y1 * si
  ypr <-  - xp.x1 * si + yp.y1 * co
  d3 <- ypr^2
  # test
  lenf <- matrix(leng, nrow=np, ncol=nl, byrow=T)
  outside <- (xpr < 0 | xpr > lenf)
  if(any(outside))
    d3[outside] <- Inf

  dmin <- matrix(pmin(d1, d2, d3),nrow=np, ncol=nl)
  return(sqrt(dmin))
}



