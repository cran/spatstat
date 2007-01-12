#
#
#    poisson.S
#
#    $Revision: 1.6 $	$Date: 2007/01/11 03:36:02 $
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
         creator  = "Poisson",
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
         irange = function(...) { 0 },
         version=versionstring.spatstat()
  )
  class(out) <- "interact"
  return(out)
}
