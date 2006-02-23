#
#   unique.ppp.R
#
# $Revision: 1.3 $  $Date: 2006/02/22 11:18:44 $
#

unique.ppp <- function(x, ...) {
  verifyclass(x, "ppp")
  dupe <- duplicated.ppp(x)
  return(x[!dupe])
}

duplicated.ppp <- function(x, ...) {
  verifyclass(x, "ppp")
  d <- pairdist(x)
  equal <- (d == 0)
  if(is.marked(x))
    equal <- equal & outer(x$marks, x$marks, "==")
  duped <- equal & (col(equal) < row(equal))
  duplic <- apply(duped, 1, any)
  return(duplic)
}

multiplicity.ppp <- function(x) {
  verifyclass(x, "ppp")
  d <- pairdist(x)
  equal <- (d == 0)
  if(is.marked(x))
    equal <- equal & outer(x$marks, x$marks, "==")
  return(as.integer(matrowsum(equal)))
}
  
  
