#
#   unique.ppp.R
#
# $Revision: 1.10 $  $Date: 2007/04/20 03:14:47 $
#

unique.ppp <- function(x, ...) {
  verifyclass(x, "ppp")
  dupe <- duplicated.ppp(x)
  return(x[!dupe])
}

duplicated.ppp <- function(x, ...) {
  verifyclass(x, "ppp")
  n <- x$n
  # marked points - split by mark value
  if(is.marked(x)) {
    m <- marks(x)
    um <- unique(m)
    xx <- unmark(x)
    result <- logical(n)
    for(i in seq(um)) {
      sub <- (m == um[i])
      result[sub] <- duplicated.ppp(xx[sub])
    }
    return(result)
  }
  # unmarked points
  # check for duplication of x and y separately (a necessary condition)
  xx <- x$x
  yy <- x$y
  possible <- duplicated(xx) & duplicated(yy)
  if(!any(possible))
    return(possible)
  # split by x coordinate of duplicated x values
  result <- possible
  xvals <- unique(xx[possible])
  for(xvalue in xvals) {
    sub <- (xx == xvalue)
    # compare y values
    result[sub] <- duplicated(yy[sub])
  }
  return(result)
}

  
multiplicity.ppp <- function(x) {
  verifyclass(x, "ppp")
  xx <- x$x
  yy <- x$y
  equal <- outer(xx, xx, "==") & outer(yy, yy, "==")
  if(is.marked(x)) {
    marx <- marks(x, dfok=FALSE)
    equal <- equal & outer(marx, marx, "==")
  }
  return(matrowsum(equal))
}
  
  
