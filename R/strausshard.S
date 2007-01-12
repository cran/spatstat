#
#
#    strausshard.S
#
#    $Revision: 2.8 $	$Date: 2007/01/11 03:36:02 $
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
         creator = "StraussHard",
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
         print = NULL,         # default OK
         interpret =  function(coeffs, self) {
           loggamma <- as.numeric(coeffs[1])
           gamma <- exp(loggamma)
           return(list(param=list(gamma=gamma),
                       inames="interaction parameter gamma",
                       printable=round(gamma,4)))
         },
         valid = function(coeffs, self) {
           gamma <- (self$interpret)(coeffs, self)$param$gamma
           return(is.finite(gamma))
         },
         project = function(coeffs, self) {
           gamma <- (self$interpret)(coeffs, self)$param$gamma
           if(is.na(gamma)) 
             coeffs[1] <- 0
           else if(!is.finite(gamma)) 
             coeffs[1] <- log(.Machine$double.xmax)
           return(coeffs)
         },
         irange = function(self, coeffs=NA, epsilon=0, ...) {
           r <- self$par$r
           hc <- self$par$hc
           if(any(is.na(coeffs)))
             return(r)
           loggamma <- coeffs[1]
           if(abs(loggamma) <= epsilon)
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
