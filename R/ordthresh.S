#
#
#    ordthresh.S
#
#    $Revision: 1.8 $	$Date: 2008/04/08 10:06:23 $
#
#    Ord process with threshold potential
#
#    OrdThresh()  create an instance of the Ord process
#                 [an object of class 'interact']
#                 with threshold potential
#	
#
# -------------------------------------------------------------------
#	

OrdThresh <- function(r) {
  out <- 
  list(
         name     = "Ord process with threshold potential",
         creator  = "OrdThresh",
         family    = ord.family,
         pot      = function(d, par) {
                         (d <= par$r)
                    },
         par      = list(r = r),
         parnames = "threshold distance",
         init     = function(self) {
                      r <- self$par$r
                      if(!is.numeric(r) || length(r) != 1 || r <= 0)
                       stop("threshold distance r must be a positive number")
                    },
         update = NULL,  # default OK
         print = NULL,    # default OK
         interpret =  function(coeffs, self) {
           loggamma <- as.numeric(coeffs[1])
           gamma <- exp(loggamma)
           return(list(param=list(gamma=gamma),
                       inames="interaction parameter gamma",
                       printable=round(gamma,4)))
         },
         irange = function(...) {
           return(Inf)
         },
       version=versionstring.spatstat()
  )
  class(out) <- "interact"
  out$init(out)
  return(out)
}
