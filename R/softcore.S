#
#
#    softcore.S
#
#    $Revision: 2.7 $   $Date: 2007/01/11 03:34:51 $
#
#    Soft core processes.
#
#    Softcore()    create an instance of a soft core process
#                 [an object of class 'interact']
#
#
# -------------------------------------------------------------------
#

Softcore <- function(kappa) {
  out <- 
  list(
         name     = "Soft core process",
         creator  = "Softcore",
         family   = pairwise.family,
         pot      = function(d, par) {
                        -d^(-2/par$kappa)
                    },
         par      = list(kappa = kappa),
         parnames = "Exponent kappa",
         init     = function(self) {
                      kappa <- self$par$kappa
                      if(!is.numeric(kappa) || length(kappa) != 1 ||
                         kappa <= 0 || kappa >= 1)
                       stop("Exponent kappa must be a positive \
number less than 1")
                    },
         update = NULL,  # default OK
         print = NULL,    # default OK
         interpret =  function(coeffs, self) {
           theta <- as.numeric(coeffs[1])
           sigma <- theta^(self$par$kappa/2)
           return(list(param=list(sigma=sigma),
                       inames="interaction parameter sigma",
                       printable=sigma))
         },
         valid = function(coeffs, self) {
           theta <- coeffs[1]
           return(is.finite(theta) && (theta >= 0))
         },
         project = function(coeffs, self) {
           theta <- coeffs[1]
           if(is.na(theta))
             theta <- 0
           coeffs[] <- max(0, theta)
           return(coeffs)
         },
         irange = function(self, coeffs=NA, epsilon=0, ...) {
           # distance d beyond which log(interaction factor) <= epsilon
           if(any(is.na(coeffs)) || epsilon == 0)
             return(Inf)
           theta <- as.numeric(coeffs[1])
           kappa <- self$par$kappa
           return((theta/epsilon)^(kappa/2))
         },
       version=versionstring.spatstat()
  )
  class(out) <- "interact"
  out$init(out)
  return(out)
}

