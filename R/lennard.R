#
#
#    lennard.R
#
#    $Revision: 1.11 $	$Date: 2007/01/11 03:36:02 $
#
#    Lennard-Jones potential
#
#
# -------------------------------------------------------------------
#	

LennardJones <- function() {
  out <- 
  list(
         name     = "Lennard-Jones potential",
         creator  = "LennardJones",
         family    = pairwise.family,
         pot      = function(d, par) {
                         array(c(d^{-12},-d^{-6}),dim=c(dim(d),2))
                    },
         par      = list(),
         parnames = character(),
         init     = function(...){}, # do nothing
         update = NULL, # default OK
         print = NULL,    # default OK
         interpret =  function(coeffs, self) {
           theta1 <- as.numeric(coeffs[1])
           theta2 <- as.numeric(coeffs[2])
           if(theta1 <= 0) {
             # Fitted regular parameter sigma^12 is negative
             sigma <- NA
             tau <- NA
           }
           else {
             sigma <- theta1^(1/12)
             tau <- theta2/sqrt(theta1)
           }
           return(list(param=list(sigma=sigma, tau=tau),
                       inames="interaction parameters",
                       printable=round(c(sigma=sigma,tau=tau),4)))
         },
         valid = function(coeffs, self) {
           p <- self$interpret(coeffs, self)$param
           return(!any(is.na(p)))
         },
         project = function(coeffs, self) {
           p <- self$interpret(coeffs, self)$param
           if(any(is.na(p)))
             stop("Don't know how to project Lennard-Jones models")
         },
         irange = function(self, coeffs=NA, epsilon=0, ...) {
           if(any(is.na(coeffs)) || epsilon == 0)
             return(Inf)
           theta1 <- abs(coeffs[1])
           theta2 <- abs(coeffs[2])
           return(max((theta1/epsilon)^(1/12), (theta2/epsilon)^(1/6)))
         },
       version=versionstring.spatstat()
  )
  class(out) <- "interact"
  out$init(out)
  return(out)
}
