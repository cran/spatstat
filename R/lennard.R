#
#
#    lennard.R
#
#    $Revision: 1.6 $	$Date: 2005/03/12 01:17:43 $
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
         family    = pairwise.family,
         pot      = function(d, par) {
                         array(c(d^{-12},-d^{-6}),dim=c(dim(d),2))
                    },
         par      = list(),
         parnames = character(),
         init     = function(...){}, # do nothing
         update = function(...){},  # default OK
         print = NULL,    # default OK
         interpret =  function(coeffs, self) {
           theta1 <- coeffs[["Interact.1"]]
           theta2 <- coeffs[["Interact.2"]]
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
           theta1 <- abs(coeffs[["Interact.1"]])
           theta2 <- abs(coeffs[["Interact.2"]])
           return(max((theta1/epsilon)^(1/12), (theta2/epsilon)^(1/6)))
         }
  )
  class(out) <- "interact"
  out$init(out)
  return(out)
}
