#
#
#    poisson.S
#
#    $Revision: 1.3 $	$Date: 2001/08/07 11:52:17 $
#
#    The Poisson process
#
#    Poisson()    create an object of class 'interact' describing
#                 the (null) interpoint interaction structure
#                 of the Poisson process.
#	
#
# -------------------------------------------------------------------
#	

Poisson <- function() {
  out <- 
  list(
         name     = "Poisson process",
         family   = NULL,
         pot      = NULL,
         par      = NULL,
         parnames = NULL,
         init     = function(...) { },
         update   = function(...) { },
         print    = function(self) {
           cat("Poisson process\n")
           invisible()
         }
  )
  class(out) <- "interact"
  return(out)
}
