#
#     dg.S
#
#    $Revision: 1.5 $	$Date: 2005/03/22 02:27:58 $
#
#     Diggle-Gratton pair potential
#
#
DiggleGratton <- function(delta, rho) {
  out <- 
  list(
         name     = "Diggle-Gratton process",
         family    = pairwise.family,
         pot      = function(d, par) {
                       delta <- par$delta
                       rho <- par$rho
                       above <- (d > rho)
                       inrange <- (!above) & (d > delta)
                       h <- above + inrange * (d - delta)/(rho - delta)
                       return(log(h))
                    },
         par      = list(delta=delta, rho=rho),
         parnames = list("lower limit delta", "upper limit rho"),
         init     = function(self) {
                      r <- self$par$delta
                      r <- self$par$rho
                      if(!is.numeric(delta) || length(delta) != 1)
                       stop("lower limit delta must be a single number")
                      if(!is.numeric(rho) || length(rho) != 1)
                       stop("upper limit rho must be a single number")
                      stopifnot(delta >= 0)
                      stopifnot(rho > delta)
                      stopifnot(is.finite(rho))
                    },
         update = NULL,  # default OK
         print = NULL,    # default OK
         interpret =  function(coeffs, self) {
           kappa <- coeffs[["Interaction"]]
           return(list(param=list(kappa=kappa),
                       inames="exponent kappa",
                       printable=round(kappa,4)))
         },
         valid = function(coeffs, self) {
           kappa <- ((self$interpret)(coeffs, self))$param$kappa
           return(is.finite(kappa) && (kappa >= 0))
         },
         project = function(coeffs, self) {
           kappa <- coeffs[["Interaction"]]
           coeffs[["Interaction"]] <-
             if(is.na(kappa)) 0 else max(0, kappa)
           return(coeffs)
         },
         irange = function(self, coeffs=NA, epsilon=0, ...) {
           rho <- self$par$rho
           if(all(is.na(coeffs)))
             return(rho)
           kappa <- coeffs[["Interaction"]]
           delta <- self$par$delta
           if(abs(kappa) <= epsilon)
             return(delta)
           else return(rho)
         }
  )
  class(out) <- "interact"
  out$init(out)
  return(out)
}
