#
#  is.subset.owin.R
#
#  $Revision: 1.2 $   $Date: 2004/02/09 13:33:09 $
#
#  Determine whether a window is a subset of another window
#
#  is.subset.owin()
#
is.subset.owin <- function(A, B) {
  A <- as.owin(A)
  B <- as.owin(B)

  if(B$type == "rectangle") {
    # Some cases can be resolved using convexity of B
    
    # (1) A is also a rectangle
   if(A$type == "rectangle") {
     xx <- A$xrange[c(1,2,2,1)]
     yy <- A$yrange[c(1,1,2,2)]
     ok <- inside.owin(xx, yy, B)
     return(all(ok))
   } 
    # (2) A is polygonal
    # Then A is a subset of B iff,
    # for every constituent polygon of A with positive sign,
    # the vertices are all in B
   if(A$type == "polygonal") {
     okpolygon <- function(a, B) {
       if(area.xypolygon(a) < 0) return(TRUE)
       ok <- inside.owin(a$x, a$y, B)
       return(all(ok))
     }
     ok <- unlist(lapply(A$bdry, okpolygon, B=B))
     return(all(ok))
   }
    # (3) Feeling lucky
    # Test whether the bounding box of A is a subset of B
    # Then a fortiori, A is a subset of B
   AA <- bounding.box(A)
   if(is.subset.owin(AA, B))
     return(TRUE)
   
 }
 # In all other cases, convexity cannot be invoked
 # Discretise
  a <- as.mask(A)
  xx <- as.vector(raster.x(a)[a$m])
  yy <- as.vector(raster.y(a)[a$m])
  ok <- inside.owin(xx, yy, B)
  return(all(ok))

}
