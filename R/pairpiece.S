#
#
#    pairpiece.S
#
#    $Revision: 1.5 $	$Date: 2004/01/08 04:50:02 $
#
#    A pairwise interaction process with piecewise constant potential
#
#    PairPiece()   create an instance of the process
#                 [an object of class 'interact']
#	
#
# -------------------------------------------------------------------
#	

PairPiece <- function(r) {
  out <- 
  list(
         name     = "Piecewise constant pairwise interaction process",
         family    = pairwise.family,
         pot      = function(d, par) {
                         r <- par$r
                         nr <- length(r)
                         out <- array(FALSE, dim=c(dim(d), nr))
                         out[,,1] <-  ifelse(d <= r[1], 1, 0)
                         for(i in 2:nr) 
                           out[,,i] <- ifelse((d > r[i-1]) & (d <= r[i]), 1, 0)
                         out
                    },
         par      = list(r = r),
         parnames = "interaction thresholds",
         init     = function(self) {
                      r <- self$par$r
                      if(!is.numeric(r) || length(r) < 2 || !all(r > 0))
                       stop("interaction thresholds r must be positive numbers")
                      if(!all(diff(r) > 0))
                        stop("interaction thresholds r must be strictly increasing")
                    },
         update = NULL,  # default OK
         print = NULL,    # default OK
         interpret =  function(coeffs, self) {
           r <- self$par$r
           npiece <- length(r)
           # extract coefficients
           vnames <- paste("Interact.", 1:npiece, sep="")
           thetas <- coeffs[vnames]
           gammas <- exp(thetas)
           # name them
           gn <- gammas
           names(gn) <- paste("(", c(0,r[-npiece]),",", r, "]", sep="")
           #
           return(list(param=list(gammas=gammas),
                       inames="interaction parameters gamma_i",
                       printable=round(gn,4)))
         },
         siminfo = NULL
  )
  class(out) <- "interact"
  out$init(out)
  return(out)
}
