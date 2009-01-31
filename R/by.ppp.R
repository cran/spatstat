#
#   by.ppp.R
#
#  $Revision: 1.3 $  $Date: 2008/07/22 21:39:42 $
#

by.ppp <- function(data, INDICES=marks(data), FUN, ...) {
  if(missing(INDICES))
    INDICES <- marks(data, dfok=FALSE)
  if(missing(FUN))
    stop("FUN is missing")
  y <- split(data, INDICES)
  z <- list()
  for(i in seq(along=y))
    z[[i]] <- FUN(y[[i]], ...)
  names(z) <- names(y)
  z <- as.listof(z)
  return(z)
}
