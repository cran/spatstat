#
#
#    strauss.S
#
#    $Revision: 1.6 $	$Date: 2002/08/05 14:18:50 $
#
#    The Strauss process
#
#    Strauss()    create an instance of the Strauss process
#                 [an object of class 'interact']
#	
#
# -------------------------------------------------------------------
#	

Strauss <- function(r) {
  out <- 
  list(
         name     = "Strauss process",
         family    = pairwise.family,
         pot      = function(d, par) {
                         d <= par$r
                    },
         par      = list(r = r),
         parnames = "interaction distance",
         init     = function(self) {
                      r <- self$par$r
                      if(!is.numeric(r) || length(r) != 1 || r <= 0)
                       stop("interaction distance r must be a positive number")
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
