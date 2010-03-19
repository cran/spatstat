#
#
#    hardcore.S
#
#    $Revision: 2.10 $	$Date: 2008/04/08 10:07:11 $
#
#    The Hard core process
#
#    Hardcore()     create an instance of the Hard Core process
#                      [an object of class 'interact']
#	
#
# -------------------------------------------------------------------
#	

Hardcore <- function(hc) {
  out <- 
  list(
         name   = "Hard core process",
         creator = "Hardcore",
         family  = pairwise.family,
         pot    = function(d, par) {
           v <- 0 * d
           v[ d <= par$hc ] <-  (-Inf)
           attr(v, "IsOffset") <- TRUE
           v
         },
         par    = list(hc = hc),
         parnames = "hard core distance", 
         init   = function(self) {
           hc <- self$par$hc
           if(!is.numeric(hc) || length(hc) != 1 || hc <= 0)
             stop("hard core distance hc must be a positive number")
         },
         update = NULL,       # default OK
         print = NULL,        # default OK
         interpret =  function(coeffs, self) {
           return(NULL)
         },
         valid = function(coeffs, self) {
           return(TRUE)
         },
         project = function(coeffs, self) {
           return(coeffs)
         },
         irange = function(self, coeffs=NA, epsilon=0, ...) {
           hc <- self$par$hc
           return(hc)
         },
       version=versionstring.spatstat()
  )
  class(out) <- "interact"
  (out$init)(out)
  return(out)
}
