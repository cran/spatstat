#
#
#    strausshard.S
#
#    $Revision: 1.5 $	$Date: 2001/08/10 09:13:59 $
#
#    The Strauss/hard core process
#
#    StraussHard()     create an instance of the Strauss-hardcore process
#                      [an object of class 'interact']
#	
#
# -------------------------------------------------------------------
#	

StraussHard <- function(r, hc) {
  out <- 
  list(
         name   = "Strauss - hard core process",
         family  = pairwise.family,
         pot    = function(d, par) {
           v <- ifelse(d <= par$r, 1, 0)
           v[ d <= par$hc ] <-  (-Inf)
           v
         },
         par    = list(r = r, hc = hc),
         parnames = c("interaction distance",
                      "hard core distance"), 
         init   = function(self) {
           r <- self$par$r
           hc <- self$par$hc
           if(!is.numeric(hc) || length(hc) != 1 || hc <= 0)
             stop("hard core distance hc must be a positive number")
           if(!is.numeric(r) || length(r) != 1 || r <= hc)
             stop("interaction distance r must be a number greater than hc")
         },
         update = NULL,       # default OK
         print = NULL         # default OK
  )
  class(out) <- "interact"
  (out$init)(out)
  return(out)
}
