#
#
#    pairpiece.S
#
#    $Revision: 1.2 $	$Date: 2002/05/13 12:41:10 $
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
         print = NULL    # default OK
  )
  class(out) <- "interact"
  out$init(out)
  return(out)
}
