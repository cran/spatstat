#
#
#    poisson.S
#
#    $Revision: 1.4 $	$Date: 2005/03/11 22:55:48 $
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
         },
         valid = function(...) { TRUE },
         project = function(coeffs, ...) { coeffs },
         irange = function(...) { 0 }
  )
  class(out) <- "interact"
  return(out)
}
