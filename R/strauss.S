#
#
#    strauss.S
#
#    $Revision: 2.1 $	$Date: 2004/01/27 08:06:38 $
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
         }
  )
  class(out) <- "interact"
  out$init(out)
  return(out)
}
