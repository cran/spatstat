#
# clickpoly.R
#
#
# $Revision: 1.1 $  $Date: 2007/03/09 04:11:39 $
#
#
clickpoly <- function(add=FALSE, n=NULL, ...) {
  if((!add) | dev.cur() == 1) {
    plot(0,0,type="n", xlab="", ylab="", xlim=c(0,1), ylim=c(0,1), asp=1.0,
         axes=FALSE)
    rect(0,0,1,1)
  }
  if(!is.null(n)) 
    cat(paste("click", n, "times in window\n"))
  else
    cat(paste("to add points: click left mouse button in window\n",
              "      to exit: click middle mouse button\n",
              "[The last point should NOT repeat the first point]\n"))
  arglist <- resolve.defaults(list(...), list(type="o"))
  if(!is.null(n)) arglist <- append(list(n=n), arglist)
  xy <- do.call("locator", arglist)
  if(area.xypolygon(xy) < 0)
    xy <- lapply(xy, rev)
  result <- owin(poly=xy)
  plot(result, add=TRUE)
  return(result)
}

  
