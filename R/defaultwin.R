#
#
#  defaultwin.R
#
#   $Revision: 1.4 $   $Date: 2009/07/17 05:57:56 $
#

default.expand <- function(object, m=2, epsilon=1e-6) {
  verifyclass(object, "ppm")
  # data window
  w <- as.owin(data.ppm(object))
  # interaction range of model
  rr <- reach(object, epsilon=epsilon)
  if(!is.finite(rr))
    return(NULL)
  if(!is.numeric(m) || length(m) != 1 || m < 1)
    stop("m should be a single number >= 1")
  mr <- m * rr
  if(w$type == "rectangle") 
    return(owin(w$xrange + c(-1,1) * mr, w$yrange + c(-1,1) * mr))
  else
    return(dilation.owin(w, mr))
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
    return(erosion.owin(w, rr))
}

  
