#
#
#    strauss.S
#
#    $Revision: 1.5 $	$Date: 2001/07/29 07:25:54 $
#
#    The Strauss process
#
#    Strauss()    create an instance of the Strauss process
#                 [an object of class 'interact']
#	
#
# -------------------------------------------------------------------
#	

Strauss <- function(r) {
  out <- 
  list(
         name     = "Strauss process",
         family    = pairwise.family,
         pot      = function(d, par) {
                         (d <= par$r)
                    },
         par      = list(r = r),
         parnames = "interaction distance",
         init     = function(self) {
                      r <- self$par$r
                      if(!is.numeric(r) || length(r) != 1 || r <= 0)
                       stop("interaction distance r must be a positive number")
                    },
         update = NULL,  # default OK
         print = NULL    # default OK
  )
  class(out) <- "interact"
  out$init(out)
  return(out)
}
