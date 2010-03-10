#
#
#    satpiece.S
#
#    $Revision: 1.9 $	$Date: 2008/04/11 15:17:30 $
#
#    Saturated pairwise interaction process with piecewise constant potential
#
#    SatPiece()   create an instance of the process
#                 [an object of class 'interact']
#	
#
# -------------------------------------------------------------------
#	

SatPiece <- function(r, sat) {
  out <- 
  list(
         name     = "piecewise constant Saturated pairwise interaction process",
         creator  = "SatPiece",
         family    = pairsat.family,
         pot      = function(d, par) {
                       r <- par$r
                       nr <- length(r)
                       out <- array(FALSE, dim=c(dim(d), nr))
                       out[,,1] <- (d < r[1])
                       if(nr > 1) {
                         for(i in 2:nr) 
                           out[,,i] <- (d >= r[i-1]) & (d < r[i])
                       }
                       out
                    },
         par      = list(r = r, sat=sat),
         parnames = c("interaction thresholds", "saturation parameters"),
         init     = function(self) {
                      r <- self$par$r
                      sat <- self$par$sat
                      if(!is.numeric(r) || !all(r > 0))
                        stop("interaction thresholds r must be positive numbers")
                      if(length(r) > 1 && !all(diff(r) > 0))
                        stop("interaction thresholds r must be strictly increasing")
                      if(!is.numeric(sat) || any(sat < 0))
                        stop("saturation parameters must be nonnegative numbers")
                      if(any(ceiling(sat) != floor(sat)))
                        warning("saturation parameter has a non-integer value")
                      if(length(sat) != length(r) && length(sat) != 1)
                        stop("vectors r and sat must have equal length")
                    },
         update = NULL,  # default OK
         print = NULL,    # default OK
         interpret =  function(coeffs, self) {
           r <- self$par$r
           npiece <- length(r)
           # extract coefficients
           gammas <- exp(as.numeric(coeffs))
           # name them
           gn <- gammas
           names(gn) <- paste("[", c(0,r[-npiece]),",", r, ")", sep="")
           #
           return(list(param=list(gammas=gammas),
                       inames="interaction parameters gamma_i",
                       printable=round(gn,4)))
         },
        valid = function(coeffs, self) {
           # interaction parameters gamma must be
           #   non-NA 
           #   finite, if sat > 0
           #   less than 1, if sat = Inf
           gamma <- (self$interpret)(coeffs, self)$param$gammas
           sat <- self$par$sat
           if(any(is.na(gamma)))
             return(FALSE)
           return(all((is.finite(gamma) | sat == 0)
                      & (gamma <= 1 | sat != Inf)))
        },
        project = function(coeffs, self){
          gamma <- (self$interpret)(coeffs, self)$param$gammas
          sat <- self$par$sat
          # convert gamma=NA to gamma = 1
          if(any(na <- is.na(gamma))) 
            coeffs[na] <- 0
          # convert gamma=Inf (where sat > 0) to a large finite value
          if(any(inf <- !na & (gamma == Inf) & (sat != 0)))
            coeffs[inf] <- log(.Machine$double.xmax)/self$par$sat
          # clip gamma to [0,1] where sat = Inf
          if(any(unsat <- (sat == Inf)))
            coeffs[unsat] <- pmin(0, coeffs[unsat])
          return(coeffs)
        },
        irange = function(self, coeffs=NA, epsilon=0, ...) {
          r <- self$par$r
          sat <- self$par$sat
          if(all(is.na(coeffs)))
            return(2 * max(r))
          gamma <- (self$interpret)(coeffs, self)$param$gammas
          gamma[is.na(gamma)] <- 1
          active <- (abs(log(gamma)) > epsilon) & (sat > 0)
          if(!any(active))
            return(0)
          else return(2 * max(r[active]))
        },
       version=versionstring.spatstat()
  )
  class(out) <- "interact"
  out$init(out)
  return(out)
}
