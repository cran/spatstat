#
#
#    strauss.S
#
#    $Revision: 2.5 $	$Date: 2006/10/10 04:22:48 $
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
         },
         valid = function(coeffs, self) {
           gamma <- ((self$interpret)(coeffs, self))$param$gamma
           return(is.finite(gamma) && (gamma <= 1))
         },
         project = function(coeffs, self) {
           loggamma <- coeffs[["Interaction"]]
           coeffs[["Interaction"]] <-
             if(is.na(loggamma)) 0 else min(0, loggamma)
           return(coeffs)
         },
         irange = function(self, coeffs=NA, epsilon=0, ...) {
           r <- self$par$r
           if(any(is.na(coeffs)))
             return(r)
           loggamma <- coeffs[["Interaction"]]
           if(abs(loggamma) <= epsilon)
             return(0)
           else
             return(r)
           
         }
  )
  class(out) <- "interact"
  out$init(out)
  return(out)
}
