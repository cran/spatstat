#
#
#    ordthresh.S
#
#    $Revision: 1.2 $	$Date: 2002/08/05 14:18:50 $
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
         family    = ord.family,
         pot      = function(d, par) {
                         ifelse(d <= par$r, 1, 0)
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
           loggamma <- coeffs[["Interaction"]]
           gamma <- exp(loggamma)
           return(list(param=list(gamma=gamma),
                       inames="interaction parameter gamma",
                       printable=round(gamma,4)))
         }
  )
  class(out) <- "interact"
  out$init(out)
  return(out)
}
