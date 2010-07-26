#
#
#    badgey.S
#
#    $Revision: 1.5 $	$Date: 2010/07/18 08:46:28 $
#
#    Hybrid Geyer process
#
#    BadGey()   create an instance of the process
#                 [an object of class 'interact']
#	
#
# -------------------------------------------------------------------
#	

BadGey <- function(r, sat) {
  out <- 
  list(
         name     = "hybrid Geyer process",
         creator  = "BadGey",
         family    = pairsat.family,
         pot      = function(d, par) {
                       r <- par$r
                       nr <- length(r)
                       out <- array(FALSE, dim=c(dim(d), nr))
                       for(i in 1:nr) 
                         out[,,i] <- (d <= r[i])
                       out
                    },
         par      = list(r = r, sat=sat),
         parnames = c("interaction radii", "saturation parameters"),
         init     = function(self) {
                      r <- self$par$r
                      sat <- self$par$sat
                      if(!is.numeric(r) || !all(r > 0))
                        stop("interaction radii r must be positive numbers")
                      if(length(r) > 1 && !all(diff(r) > 0))
                        stop("interaction radii r must be strictly increasing")
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
           names(gn) <- paste("[0,", r, ")", sep="")
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
       version=versionstring.spatstat(),
       # fast evaluation is available for the border correction only
       can.do.fast=function(X,correction,par) {
         return(all(correction %in% c("border", "none")))
       },
       fasteval=function(X,U,EqualPairs,pairpot,potpars,correction,
                         ..., halfway=FALSE) {
         # fast evaluator for BadGey interaction
         if(!all(correction %in% c("border", "none")))
           return(NULL)
         if(spatstat.options("fasteval") == "test")
           message("Using fast eval for BadGey")
         r  <- potpars$r
         hc <- potpars$sat
         # ensure r and sat have equal length
         if(length(r) != length(sat)) {
           if(length(r) == 1)
             r <- rep(r, length(sat))
           else if(length(sat) == 1)
             sat <- rep(sat, length(r))
           else stop("lengths or r and sat do not match")
         }
         nterms <- length(r)
         answer <- matrix(, nrow=npoints(U), ncol=nterms)
         for(k in 1:nterms) {
           # first determine saturated pair counts
           counts <- strausscounts(U, X, r[k], EqualPairs) 
           satcounts <- pmin(sat[k], counts)
           # trapdoor used by suffstat() 
           if(halfway) 
             answer[,k] <- satcounts
           else if(sat[k] == Inf)
             answer[,k] <- 2 * satcounts
           else {
             # extract counts for data points
             Uindex <- EqualPairs[,2]
             Xindex <- EqualPairs[,1]
             Xcounts <- integer(npoints(X))
             Xcounts[Xindex] <- counts[Uindex]
             # evaluate change in saturated counts of other data points
             change <- geyercounts(U, X, r[k], sat[k], Xcounts, EqualPairs)
             answer[,k] <- satcounts + change
           }
         }
         return(answer)
       }
  )
  class(out) <- "interact"
  out$init(out)
  return(out)
}
