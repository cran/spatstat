#
#    rescue.rectangle.R
# 
#    $Revision: 1.1 $   $Date: 2005/03/02 23:08:49 $
#
rescue.rectangle <- function(W) {
  verifyclass(W, "owin")

  if(W$type == "mask" && all(W$m))
     return(owin(W$xrange, W$yrange))

  if(W$type == "polygonal" && length(W$bdry) == 1) {
    x <- W$bdry[[1]]$x
    y <- W$bdry[[1]]$y
    if(length(x) == 4 && length(y) == 4 &&
       length(ux <- unique(x)) == 2 && 
       length(uy <- unique(y)) == 2)
      return(owin(ux,uy))
  }

  return(W)
}

