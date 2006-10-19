#
#   unique.ppp.R
#
# $Revision: 1.6 $  $Date: 2006/10/10 04:22:48 $
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
  ism <- is.marked(x, dfok=TRUE)
  if(ism) {
    mm <- marks(x, dfok=TRUE)
    possible <- possible & duplicated(mm)
  }
  if(!any(possible))
    return(duped)
  crossmatch <- function(a, j) { outer(a[j], a[!j], "==") }
  equal <- crossmatch(xx, possible)
  equal <- equal & crossmatch(yy, possible)
  if(ism) {
    if(!is.data.frame(mm))
      equal <- equal & crossmatch(mm, possible)
    else {
      for(i in seq(ncol(mm)))
        equal <- equal & crossmatch(mm[, i], possible)
    }
  }
  duped[possible] <- apply(equal, 1, any)
  return(duped)
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
  
  
