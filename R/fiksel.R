#
#
#    fiksel.R
#
#    $Revision: 1.1 $	$Date: 2010/05/14 08:59:29 $
#
#    Fiksel interaction 
#    
#    ee Stoyan Kendall Mcke 1987 p 161
#
# -------------------------------------------------------------------
#	

Fiksel <- function(r, hc, kappa) {
  out <- 
  list(
         name   = "Fiksel process",
         creator = "Fiksel",
         family  = pairwise.family,
         pot    = function(d, par) {
           v <- (d <= par$r) * exp( - d * par$kappa)
           v[ d <= par$hc ] <-  (-Inf)
           v
         },
         par    = list(r = r, hc = hc, kappa=kappa),
         parnames = c("interaction distance",
                      "hard core distance",
                      "rate parameter"), 
         init   = function(self) {
           r <- self$par$r
           hc <- self$par$hc
           kappa <- self$par$kappa
           if(!is.numeric(hc) || length(hc) != 1 || hc <= 0)
             stop("hard core distance hc must be a positive number")
           if(!is.numeric(r) || length(r) != 1 || r <= hc)
             stop("interaction distance r must be a number greater than hardcore dstance hc")
           if(!is.numeric(kappa) || length(kappa) != 1)
             stop("rate parameter kappa must be a single number")
         },
         update = NULL,       # default OK
         print = NULL,         # default OK
         interpret =  function(coeffs, self) {
           a <- as.numeric(coeffs[1])
           return(list(param=list(a=a),
                       inames="interaction strength a",
                       printable=round(a,2)))
         },
         valid = function(coeffs, self) {
           a <- (self$interpret)(coeffs, self)$param$a
           return(is.finite(a))
         },
         project = function(coeffs, self) {
           a <- (self$interpret)(coeffs, self)$param$a
           if(is.na(a)) 
             coeffs[1] <- 0
           else if(!is.finite(a)) 
             coeffs[1] <- log(.Machine$double.xmax)
           return(coeffs)
         },
         irange = function(self, coeffs=NA, epsilon=0, ...) {
           r <- self$par$r
           hc <- self$par$hc
           if(any(is.na(coeffs)))
             return(r)
           a <- coeffs[1]
           if(abs(a) <= epsilon)
             return(hc)
           else
             return(r)
         },
       version=versionstring.spatstat()
  )
  class(out) <- "interact"
  (out$init)(out)
  return(out)
}
