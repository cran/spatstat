#
#
#    strausshard.S
#
#    $Revision: 1.7 $	$Date: 2003/03/12 05:24:28 $
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
         print = NULL,         # default OK
         interpret =  function(coeffs, self) {
           loggamma <- coeffs[["Interaction"]]
           gamma <- exp(loggamma)
           return(list(param=list(gamma=gamma),
                       inames="interaction parameter gamma",
                       printable=round(gamma,4)))
         },
         rmhmodel = function(fit, self) {
           # not implemented if there is a trend 
           if(!no.trend.ppm(fit))
             return(NULL)
           coeffs <- fit$coef
           beta <- exp(coeffs[["(Intercept)"]])
           gamma <- exp(coeffs[["Interaction"]])
           r <- self$par$r
           hc <- self$par$hc
           w <- fit$Q$data$window
           return(list(cif='straush',
                       par=c(beta=beta,gamma=gamma,r=r,hc=hc),
                       w=w))
         }
  )
  class(out) <- "interact"
  (out$init)(out)
  return(out)
}
