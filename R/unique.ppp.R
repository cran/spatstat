#
#   unique.ppp.R
#
# $Revision: 1.5 $  $Date: 2006/03/28 10:03:43 $
#

unique.ppp <- function(x, ...) {
  verifyclass(x, "ppp")
  dupe <- duplicated.ppp(x)
  return(x[!dupe])
}

duplicated.ppp <- function(x, ...) {
  verifyclass(x, "ppp")
  xx <- x$x
  yy <- x$y
  duped <- rep(FALSE, x$n)
  possible <- duplicated(xx) & duplicated(yy)
  if(is.marked(x)) 
    possible <- possible & duplicated((mm <- x$marks))
  if(!any(possible))
    return(duped)
  equal <- outer(xx[possible], xx[!possible], "==")
  equal <- equal & outer(yy[possible], yy[!possible], "==")
  if(is.marked(x))
    equal <- equal & outer(mm[possible], mm[!possible], "==")
  duped[possible] <- apply(equal, 1, any)
  return(duped)
}

multiplicity.ppp <- function(x) {
  verifyclass(x, "ppp")
  xx <- x$x
  yy <- x$y
  equal <- outer(xx, xx, "==") & outer(yy, yy, "==")
  if(is.marked(x))
    equal <- equal & outer(x$marks, x$marks, "==")
  return(matrowsum(equal))
}
  
  
