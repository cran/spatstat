#
#    rescue.rectangle.R
# 
#    $Revision: 1.5 $   $Date: 2007/10/24 09:41:15 $
#
rescue.rectangle <- function(W) {
  verifyclass(W, "owin")

  if(W$type == "mask" && all(W$m))
     return(owin(W$xrange, W$yrange, unitname=unitname(W)))

  if(W$type == "polygonal" && length(W$bdry) == 1) {
    x <- W$bdry[[1]]$x
    y <- W$bdry[[1]]$y
    if(length(x) == 4 && length(y) == 4) {
      # could be a rectangle
      veryunique <- function(z) {
        uz <- sort(unique(z))
        close <- (diff(uz) < 2 * .Machine$double.eps)
        uz <- uz[c(TRUE, !close)]
        return(uz)
      }
      ux <- veryunique(x)
      uy <- veryunique(y)
      if(length(ux) == 2 && length(uy) == 2)
        return(owin(ux,uy, unitname=unitname(W)))
    }
  }
  
  return(W)
}

