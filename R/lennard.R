#
#
#    lennard.R
#
#    $Revision: 1.13 $	$Date: 2010/06/04 04:31:23 $
#
#    Lennard-Jones potential
#
#
# -------------------------------------------------------------------
#	

LennardJones <- function() {
  out <- 
  list(
         name     = "Lennard-Jones process",
         creator  = "LennardJones",
         family    = pairwise.family,
         pot      = function(d, par) {
                         d6 <- d^{-6}
                         p <- array(c(-d6^2,d6),dim=c(dim(d),2))
                         return(p)
                    },
         par      = list(),
         parnames = character(),
         init     = function(...){}, # do nothing
         update = NULL, # default OK
         print = NULL,    # default OK
         interpret =  function(coeffs, self) {
           theta1 <- as.numeric(coeffs[1])
           theta2 <- as.numeric(coeffs[2])
           if(sign(theta1) * sign(theta2) == 1) {
             sigma <- (theta1/theta2)^(1/6)
             epsilon <- (theta2^2)/(4 * theta1)
           } else {
             sigma <- NA
             epsilon <- NA
           }
           return(list(param=list(sigma=sigma, epsilon=epsilon),
                       inames="interaction parameters",
                       printable=round(c(sigma=sigma,epsilon=epsilon),4)))
         },
         valid = function(coeffs, self) {
           p <- self$interpret(coeffs, self)$param
           return(all(!is.na(p) & (p > 0)))
         },
         project = function(coeffs, self) {
           if(!self$valid(coeffs, self))
             coeffs[] <- 0
           return(coeffs)
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
