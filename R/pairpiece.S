#
#
#    pairpiece.S
#
#    $Revision: 1.9 $	$Date: 2005/03/12 01:10:07 $
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
                       out[,,1] <-  ifelse(d < r[1], 1, 0)
                       if(nr > 1) {
                         for(i in 2:nr) 
                           out[,,i] <- ifelse((d >= r[i-1]) & (d < r[i]), 1, 0)
                       }
                       out
                    },
         par      = list(r = r),
         parnames = "interaction thresholds",
         init     = function(self) {
                      r <- self$par$r
                      if(!is.numeric(r) || !all(r > 0))
                       stop("interaction thresholds r must be positive numbers")
                      if(length(r) > 1 && !all(diff(r) > 0))
                        stop("interaction thresholds r must be strictly increasing")
                    },
         update = NULL,  # default OK
         print = NULL,    # default OK
         interpret =  function(coeffs, self) {
           r <- self$par$r
           npiece <- length(r)
           # extract coefficients
           vnames <- if(npiece == 1) "Interaction" else
                      paste("Interact.", 1:npiece, sep="")
           thetas <- coeffs[vnames]
           gammas <- exp(thetas)
           # name them
           gn <- gammas
           names(gn) <- paste("[", c(0,r[-npiece]),",", r, ")", sep="")
           #
           return(list(param=list(gammas=gammas),
                       inames="interaction parameters gamma_i",
                       printable=round(gn,4)))
         },
        valid = function(coeffs, self) {
           # interaction parameters gamma
           gamma <- (self$interpret)(coeffs, self)$param$gammas
           return(all(gamma <= 1) || gamma[1] == 0)
        },
        project = function(coeffs, self){
           # interaction parameters gamma
           gamma <- (self$interpret)(coeffs, self)$param$gammas
           if(all(gamma <= 1))
             return(coeffs)
           # clip to 1
           r <- self$par$r
           npiece <- length(r)
           vnames <- if(npiece == 1) "Interaction" else
                     paste("Interact.", 1:npiece, sep="")
           coeffs[vnames] <- pmin(0, coeffs[vnames])
           return(coeffs)
        },
        irange = function(self, coeffs=NA, epsilon=0, ...) {
          r <- self$par$r
          if(is.na(coeffs))
            return(max(r))
          gamma <- (self$interpret)(coeffs, self)$param$gammas
          gamma[is.na(gamma)] <- 1
          active <- (abs(log(gamma)) > epsilon)
          if(!any(active))
            return(0)
          else return(max(r[active]))
        }
  )
  class(out) <- "interact"
  out$init(out)
  return(out)
}
