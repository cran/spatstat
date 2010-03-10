#
#   by.ppp.R
#
#  $Revision: 1.4 $  $Date: 2010/03/08 08:23:04 $
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
