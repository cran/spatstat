#
#
#    geyer.S
#
#    $Revision: 1.3 $	$Date: 2002/08/05 14:18:50 $
#
#    Geyer's saturation process
#
#    Geyer()    create an instance of Geyer's saturation process
#                 [an object of class 'interact']
#
# ------------------------------------------------------------------
#    Note: if you want to imitate this, remember that 'pairsat.family'
#    expects the saturation parameter 'sat' to be called $par$saturate
#    in this 'interact' object.
# -------------------------------------------------------------------
#	

Geyer <- function(r, sat) {
  out <- 
  list(
         name     = "Geyer saturation process",
         family    = pairsat.family,
         pot      = function(d, par) {
                         ifelse(d <= par$r, 1, 0)  # same as for Strauss
                    },
         par      = list(r = r, saturate=sat),
         parnames = c("interaction distance","saturation parameter"),
         init     = function(self) {
                      r <- self$par$r
                      sat <- self$par$saturate
                      if(!is.numeric(r) || length(r) != 1 || r <= 0)
                       stop("interaction distance r must be a positive number")
                      if(!is.numeric(sat) || length(sat) != 1 || sat < 1)
                       stop("saturation parameter sat must be a number >= 1")
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
