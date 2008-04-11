#
#
#    geyer.S
#
#    $Revision: 2.13 $	$Date: 2008/04/08 10:06:39 $
#
#    Geyer's saturation process
#
#    Geyer()    create an instance of Geyer's saturation process
#                 [an object of class 'interact']
#
#	

Geyer <- function(r, sat) {
  out <- 
  list(
         name     = "Geyer saturation process",
         creator  = "Geyer",
         family    = pairsat.family,
         pot      = function(d, par) {
                         (d <= par$r)  # same as for Strauss
                    },
         par      = list(r = r, sat=sat),
         parnames = c("interaction distance","saturation parameter"),
         init     = function(self) {
                      r <- self$par$r
                      sat <- self$par$sat
                      if(!is.numeric(r) || length(r) != 1 || r <= 0)
                       stop("interaction distance r must be a positive number")
                      if(!is.numeric(sat) || length(sat) != 1 || sat < 1)
                       stop("saturation parameter sat must be a number >= 1")
                      if(ceiling(sat) != floor(sat))
                        warning(paste("saturation parameter sat",
                                      "has a non-integer value"))
                    },
         update = NULL, # default OK
         print = NULL,    # default OK
         interpret =  function(coeffs, self) {
           loggamma <- as.numeric(coeffs[1])
           gamma <- exp(loggamma)
           return(list(param=list(gamma=gamma),
                       inames="interaction parameter gamma",
                       printable=round(gamma,4)))
         },
         valid = function(coeffs, self) {
           gamma <- (self$interpret)(coeffs, self)$param$gamma
           return(is.finite(gamma))
         },
         project = function(coeffs, self) {
           gamma <- (self$interpret)(coeffs, self)$param$gamma
           if(is.na(gamma)) 
             coeffs[1] <- 0
           else if(!is.finite(gamma)) 
             coeffs[1] <- log(.Machine$double.xmax)/self$par$sat
           return(coeffs)
         },
         irange = function(self, coeffs=NA, epsilon=0, ...) {
           r <- self$par$r
           if(any(!is.na(coeffs))) {
             loggamma <- coeffs[1]
             if(!is.na(loggamma) && (abs(loggamma) <= epsilon))
               return(0)
           }
           return(2 * r)
         },
       version=versionstring.spatstat()
  )
  class(out) <- "interact"
  out$init(out)
  return(out)
}
