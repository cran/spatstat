#
#
#    softcore.S
#
#    $Revision: 2.11 $   $Date: 2012/08/22 11:53:54 $
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
         sig0 <- par$sigma0
         if(is.na(sig0)) {
           p <- -d^(-2/par$kappa)
         } else {
           # expand around sigma0 and set large numbers to Inf
           drat <- d/sig0
           p <- -drat^(-2/par$kappa)
           p[p < -25] <- -Inf
         }
         return(p)
       },
       par      = list(kappa = NULL, sigma0=NA),  # filled in later
       parnames = c("Exponent kappa", "Initial approximation to sigma"),
       selfstart = function(X, self) {
         # self starter for Softcore
         kappa <- self$par$kappa
         # attempt to set value of 'sigma0'
         if(!is.na(self$par$sigma0)) {
           # value fixed by user or previous invocation
           return(self)
         }
         if(npoints(X) < 2) {
           # not enough points
           return(self)
         }
         s0 <- min(nndist(X))
         if(s0 == 0) {
           warning(paste("Pattern contains duplicated points:",
                         "impossible under Softcore model"))
           s0 <- mean(nndist(X))
           if(s0 == 0)
             return(self)
         }
         Softcore(kappa=kappa, sigma0=s0)           
       },
       init     = function(self) {
         kappa <- self$par$kappa
         if(!is.numeric(kappa) || length(kappa) != 1 ||
            kappa <= 0 || kappa >= 1)
           stop(paste("Exponent kappa must be a",
                      "positive number less than 1"))
       },
       update = NULL,  # default OK
       print = NULL,    # default OK
       interpret =  function(coeffs, self) {
         theta <- as.numeric(coeffs[1])
         sigma <- theta^(self$par$kappa/2)
         if(!is.na(sig0 <- self$par$sigma0))
           sigma <- sigma * sig0
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
         sig0  <- self$par$sigma0
         if(is.na(sig0)) sig0 <- 1
         return(sig0 * (theta/epsilon)^(kappa/2))
       },
       Mayer=function(coeffs, self) {
         # second Mayer cluster integral
         kappa <- self$par$kappa
         sigma <- (self$interpret)(coeffs, self)$param$sigma
         return(pi * (sigma^2) * gamma(1 - kappa))
       },
       version=NULL # filled in later
  )
  class(BlankSoftcore) <- "interact"

  Softcore <- function(kappa, sigma0=NA) {
    instantiate.interact(BlankSoftcore, list(kappa=kappa, sigma0=sigma0))
  }

  Softcore
})

                  
