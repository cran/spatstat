#
#
#    softcore.S
#
#    $Revision: 2.9 $   $Date: 2012/01/18 10:44:55 $
#
#    Soft core processes.
#
#    Softcore()    create an instance of a soft core process
#                 [an object of class 'interact']
#
#
# -------------------------------------------------------------------
#

Softcore <- local({

  BlankSoftcore <- 
  list(
         name     = "Soft core process",
         creator  = "Softcore",
         family   = "pairwise.family",  # evaluated later
         pot      = function(d, par) {
                        -d^(-2/par$kappa)
                    },
         par      = list(kappa = NULL),  # filled in later
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
           if((self$valid)(coeffs, self)) return(NULL) else return(Poisson())
         },
         irange = function(self, coeffs=NA, epsilon=0, ...) {
           # distance d beyond which log(interaction factor) <= epsilon
           if(any(is.na(coeffs)) || epsilon == 0)
             return(Inf)
           theta <- as.numeric(coeffs[1])
           kappa <- self$par$kappa
           return((theta/epsilon)^(kappa/2))
         },
       version=NULL # filled in later
  )
  class(BlankSoftcore) <- "interact"

  Softcore <- function(kappa) {
    instantiate.interact(BlankSoftcore, list(kappa=kappa))
  }

  Softcore
})

                  
