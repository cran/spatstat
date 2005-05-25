#
#
#  defaultwin.R
#
#   $Revision: 1.1 $   $Date: 2005/05/11 19:49:02 $
#

default.expand <- function(object, m=2, epsilon=1e-6) {
  verifyclass(object, "ppm")
  # data window
  w <- as.owin(data.ppm(object))
  # interaction range of model
  rr <- reach(object, epsilon=epsilon)
  if(!is.finite(rr))
    return(NULL)
  return(owin(w$xrange + c(-1,1) * m * rr,
              w$yrange + c(-1,1) * m * rr))
}

default.clipwindow <- function(object, epsilon=1e-6) {
  verifyclass(object, "ppm")
  # data window
  w <- as.owin(data.ppm(object))
  # interaction range of model
  rr <- reach(object, epsilon=epsilon)
  if(!is.finite(rr))
    return(NULL)
  if(rr == 0)
    return(w)
  else
    return(erode.owin(w, rr))
}

  
