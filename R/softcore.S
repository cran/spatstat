#
#
#    softcore.S
#
#    $Revision: 1.2 $   $Date: 2000/07/11 10:52:01 $
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
         print = NULL    # default OK
  )
  class(out) <- "interact"
  out$init(out)
  return(out)
}

