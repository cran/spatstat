#
#
#    lennard.R
#
#    $Revision: 1.2 $	$Date: 2004/01/08 04:52:40 $
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
             warning("Fitted regular parameter sigma^12 is negative")
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
         siminfo = NULL
  )
  class(out) <- "interact"
  out$init(out)
  return(out)
}
