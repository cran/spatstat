#
#
#    strauss.S
#
#    $Revision: 1.8 $	$Date: 2003/03/12 05:16:45 $
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
         rmhmodel = function(fit, self) {
           # not implemented if there is a trend 
           if(!no.trend.ppm(fit))
             return(NULL)
           coeffs <- fit$coef
           beta <- exp(coeffs[["(Intercept)"]])
           gamma <- exp(coeffs[["Interaction"]])
           r <- self$par$r
           w <- fit$Q$data$window
           return(list(cif='strauss',
                       par=c(beta=beta,gamma=gamma,r=r),
                       w=w))
         }
  )
  class(out) <- "interact"
  out$init(out)
  return(out)
}
