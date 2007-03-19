#
#
#   rmhmodel.R
#
#   $Revision: 1.21 $  $Date: 2007/03/13 03:41:56 $
#
#

rmhmodel <- function(...) {
  UseMethod("rmhmodel")
}

rmhmodel.rmhmodel <- function(model, ...) {
  return(model)
}

rmhmodel.list <- function(model, ...) {
  argnames <- c("cif","par","w","trend","types")
  ok <- argnames %in% names(model)
  do.call("rmhmodel.default", model[argnames[ok]])
}

rmhmodel.default <- function(...,
                     cif=NULL, par=NULL, w=NULL, trend=NULL, types=NULL)
{
  if(length(list(...)) > 0)
    stop(paste("Syntax should be rmhmodel(cif, par, w, trend, types)\n",
               "with arguments given by name if they are present"))
  
  # Validate parameters
  if(is.null(cif)) stop("cif is missing")
  if(is.null(par)) stop("par is missing")

  if(!is.null(w))
    w <- as.owin(w)
  
  if(!is.character(cif) || length(cif) != 1)
    stop("cif should be a character string")

  # Check that this is a recognised model
  # and look up the rules for this model
  rules <- .Spatstat.RmhTable[[cif]]
  if(is.null(rules))
    stop(paste("Unrecognised cif:", sQuote(cif)))
  
  # Turn the name of the cif into a number (for use in Fortran)
  fortran.id <- rules$fortran.id

  # Validate the model parameters and reformat them
  check <- rules$parhandler
  fortran.par <-
    if(!rules$multitype)
      check(par)
    else if(!is.null(types))
      check(par, types)
    else 
      # types vector not given - defer checking
      NULL

  # ensure it's a numeric vector
  fortran.par <- unlist(fortran.par)

  # Check that cif executable is actually loaded
  if(!is.na(fortran.id) && !is.loaded(cif))
    stop(paste("executable is not loaded for cif", sQuote(cif)))

  # Calculate reach of model
  mreach <- rules$reach(par)

###################################################################
# return augmented list  
  out <- list(cif=cif,
              par=par,
              w=w,
              trend=trend,
              types=types,
              fortran.id=fortran.id,
              fortran.par=fortran.par,
              check= if(is.null(fortran.par)) check else NULL,
              multitype.interact=rules$multitype,
              need.aux=rules$need.aux,
              reach=mreach
              )
  class(out) <- c("rmhmodel", class(out))
  return(out)
  }

print.rmhmodel <- function(x, ...) {
  verifyclass(x, "rmhmodel")

  cat("Metropolis-Hastings algorithm, model parameters\n")

  cat(paste("Conditional intensity: cif=", x$cif, "\n"))

  if(!is.null(x$types)) {
    if(length(x$types) == 1)
      cat("Univariate process.\n")
    else {
      cat("Multitype process with types =\n")
      print(x$types)
      if(!x$multitype.interact)
        cat("Interaction does not depend on type\n")
    }
  } else if(x$multitype.interact) 
    cat("Multitype process, types not yet specified.\n")
  
  cat("Numerical parameters: par =\n")
  print(x$par)
  if(is.null(x$fortran.par))
    cat("Parameters have not yet been checked for compatibility with types.\n")
  if(is.owin(x$w)) print(x$w) else cat("Window: not specified.\n")
  cat("Trend: ")
  if(!is.null(x$trend)) print(x$trend) else cat("none.\n")

}

reach.rmhmodel <- function(x, ...) {
  if(length(list(...)) == 0)
    return(x$reach)
  rules <- .Spatstat.RmhTable[[x$cif]]
  return(rules$reach(x$par, ...))
}

#####  Table of rules for handling rmh models ##################

.Spatstat.RmhTable <-
  list(
#
# 0. Poisson (special case)
#
       'poisson'=
       list(
            fortran.id=NA,
            multitype=FALSE,
            need.aux=FALSE,
            parhandler=function(par, ...) {
              pnames <- names(par)
              ok <- ("beta" %in% pnames)
              if(!any(ok))
                stop("For the Poisson process, must specify beta")
              if(any(!ok))
                stop(paste("Unrecognised",
                           ngettext(sum(!ok), "parameter", "parameters"),
                           paste(sQuote(pnames[!ok]), collapse=", "),
                           "for Poisson process"))
              if(par[["beta"]] < 0)
                stop("Negative value of beta for Poisson process")
              return(par)
            },
            reach = function(par, ...) { return(0) }
            ),
#       
# 1. Strauss.
#       
       'strauss'=
       list(
            fortran.id=1,
            multitype=FALSE,
            need.aux=FALSE,
            parhandler=function(par, ...) {
              nms <- c("beta","gamma","r")
              if(sum(!is.na(match(names(par),nms))) != 3) {
		stop(paste("For the strauss cif, par must be a named vector\n",
                           "with components beta, gamma, and r.\n",
                           "Bailing out."))
              }
              if(any(par<0))
		stop("Negative parameters for strauss cif.")
              if(par["gamma"] > 1)
		stop("For Strauss processes, gamma must be <= 1.")
              par <- par[nms]
# Square the radius to avoid squaring in the Fortran code.
              par[[3]] <- par[[3]]^2
              return(par)
            },
            reach = function(par, ...) {
              r <- par[["r"]]
              g <- par[["gamma"]]
              return(if(g == 1) 0 else r)
            }
            ),
#       
# 2. Strauss with hardcore.
#       
       'straush' =
       list(
            fortran.id=2,
            multitype=FALSE,
            need.aux=FALSE,
            parhandler=function(par, ...) {
              nms <- c("beta","gamma","r","hc")
              if(sum(!is.na(match(names(par),nms))) != 4) {
		stop(paste("For the straush cif, par must be a named vector\n",
                           "with components beta, gamma, r, and hc.\n",
                           "Bailing out."))
              }
              if(any(par<0))
		stop("Negative parameters for straush cif.")
              par <- par[nms]
# Square the radii to avoid squaring in the Fortran code.
	      par[[3]] <- par[[3]]^2
	      par[[4]] <- par[[4]]^2
              return(par)
            },
            reach = function(par, ...) {
              h <- par[["hc"]]
              r <- par[["r"]]
              g <- par[["gamma"]]
              return(if(g == 1) h else r)
            }
            ),
#       
# 3. Softcore.
#
       'sftcr' =
       list(
            fortran.id=3,
            multitype=FALSE,
            need.aux=FALSE,
            parhandler=function(par, ...) {
              	nms <- c("beta","sigma","kappa")
                if(sum(!is.na(match(names(par),nms))) != 3) {
                  stop(paste("For the sftcr cif, par must be a named vector\n",
                             "with components beta, sigma, and kappa.\n",
                             "Bailing out."))
                }
                if(any(par<0))
                  stop("Negative  parameters for sftcr cif.")
                if(par["kappa"] > 1)
                  stop("For Softcore processes, kappa must be <= 1.")
                par <- par[nms]
                return(par)
            },
            reach = function(par, ..., epsilon=0) {
              if(epsilon==0)
                return(Inf)
              kappa <- par[["kappa"]]
              sigma <- par[["sigma"]]
              return(sigma/(epsilon^(kappa/2)))
            }                        
            ),
#       
# 4. Marked Strauss.
#       
       'straussm' =
       list(
            fortran.id=4,
            multitype=TRUE,
            need.aux=FALSE,
            parhandler=function(par, types) {
              nms <- c("beta","gamma","radii")
              if(!is.list(par) || sum(!is.na(match(names(par),nms))) != 3) 
		stop(paste("For the straussm cif, par must be a named list\n",
                           "with components beta, gamma, and radii.\n",
                           "Bailing out."))
              beta <- par$beta
              if(any(is.na(beta)))
		stop("Missing values not allowed in beta.")
              ntypes <- length(types)
              if(length(beta) != ntypes)
		stop("Length of beta does not match ntypes.")
              gamma <- par$gamma
              if(!is.matrix(gamma) || sum(dim(gamma) == ntypes) != 2)
		stop("Component gamma of par is of wrong shape.")
	      gamma[is.na(gamma)] <- 0
              r <- par$radii
              if(!is.matrix(r) || sum(dim(r) == ntypes) != 2)
		stop("Component r of par is of wrong shape.")
	      r[is.na(r)] <- 0
              gamma <- t(gamma)[row(gamma)>=col(gamma)]
              r <- t(r)[row(r)>=col(r)]
# Pack up just to check for negative values:
              par <- c(beta,gamma,r)
              if(any(par<0))
		stop("Negative  parameters for straussm cif.")
# Square the radii to avoid squaring in Fortran code:
	      r <- r^2
# Repack
	      par <- c(beta,gamma,r)
              return(par)
            }, 
            reach = function(par, ...) {
              r <- par$radii
              g <- par$gamma
              operative <- ! (is.na(r) | (g == 1))
              return(max(0, r[operative]))
            }
            ),
#       
# 5. Marked Strauss with hardcore.
#       
       'straushm' = 
       list(
            fortran.id=5,
            multitype=TRUE,
            need.aux=FALSE,
            parhandler=function(par, types) {
              nms <- c("beta","gamma","iradii","hradii")
              if(!is.list(par) || sum(!is.na(match(names(par),nms))) != 4) {
		stop(paste("For the straushm cif, par must be a named list",
                           "with components beta, gamma, iradii, and hradii.",
                           "\nBailing out."))
              }
              beta <- par$beta
              ntypes <- length(types)
              if(length(beta) != ntypes)
		stop("Length of beta does not match ntypes.")
              if(any(is.na(beta)))
		stop("Missing values not allowed in beta.")
              
              gamma <- par$gamma
              if(!is.matrix(gamma) || sum(dim(gamma) == ntypes) != 2)
		stop("Component gamma of par is of wrong shape.")
              gamma[is.na(gamma)] <- 1

              iradii <- par$iradii
              if(!is.matrix(iradii) || sum(dim(iradii) == ntypes) != 2)
		stop("Component iradii of par is of wrong shape.")
              iradii[is.na(iradii)] <- 0

              hradii <- par$hradii
              if(!is.matrix(hradii) || sum(dim(hradii) == ntypes) != 2)
		stop("Component hradii of par is of wrong shape.")
              hradii[is.na(hradii)] <- 0

              gamma <- t(gamma)[row(gamma)>=col(gamma)]
              iradii <- t(iradii)[row(iradii)>=col(iradii)]
              hradii <- t(hradii)[row(hradii)>=col(hradii)]

# Pack up just to check for negative values.
              par <- c(beta,gamma,iradii,hradii)
              if(any(par<0))
                  stop("Some parameters negative.")
# Square the radii to avoid squaring them in the Fortran code.
	      iradii <- iradii^2
	      hradii <- hradii^2
# Repack.
              par <- c(beta,gamma,iradii,hradii)
              return(par)
            },
            reach=function(par, ...) {
              r <- par$iradii
              h <- par$hradii
              g <- par$gamma
              roperative <- ! (is.na(r) | (g == 1))
              hoperative <- ! is.na(h)
              return(max(0, r[roperative], h[hoperative]))
            }
            ),
#       
# 6. Diggle-Gates-Stibbard interaction
#    (function number 1 from Diggle, Gates, and Stibbard)
       
       'dgs' = 
       list(
            fortran.id=6,
            multitype=FALSE,
            need.aux=FALSE,
            parhandler=function(par, ...) {
              nms <- c("beta","rho")
              if(sum(!is.na(match(names(par),nms))) != 2) {
		stop(paste("For the dgs cif, par must be a named vector\n",
                           "with components beta and rho.\n",
                           "Bailing out."))
              }
              if(any(par<0))
		stop("Negative parameters for dgs cif.")
              par <- par[nms]
# Now tack on rho-squared to avoid squaring in the Fortran code.
	      par <- c(par,par[[2]]^2)
              return(par)
            },
            reach=function(par, ...) {
              return(par[["rho"]])
            }
            ),
#
# 7. Diggle-Gratton interaction 
#    (function number 2 from Diggle, Gates, and Stibbard).

       'diggra' =
       list(
            fortran.id=7,
            multitype=FALSE,
            need.aux=FALSE,
            parhandler=function(par, ...) {
              nms <- c("beta","kappa","delta","rho")
              if(sum(!is.na(match(names(par),nms))) != 4) {
		stop(paste("For the diggra cif, par must be a named vector\n",
                           "with components beta, kappa, delta, and rho.\n",
                           "Bailing out."))
              }
              if(any(par<0))
		stop("Negative parameters for diggra cif.")
              if(par["delta"] >= par["rho"])
		stop("Radius delta must be less than radius rho.")
              par <- par[nms]
# Now tack on delta-squared, rho-squared, and log(rho-delta)
# to avoid calculating them in  the Fortran code.
	      par <- c(par,par[[3]]^2,par[[4]]^2,log(par[[4]]-par[[3]]))
              return(par)
            },
            reach=function(par, ...) {
              return(par[["rho"]])
            }
            ),
#       
# 8. Geyer saturation model
#       
       'geyer' = 
       list(
            fortran.id=8,
            multitype=FALSE,
            need.aux=TRUE,
            parhandler=function(par, ...) {
              nms <- c("beta","gamma","r","sat")
              if(sum(!is.na(match(names(par),nms))) != 4) {
		stop(paste("For the geyer cif, par must be a named vector\n",
                           "with components beta, gamma, r, and sat.\n",
                           "Bailing out."))
              }
              if(any(par<0))
		stop("Negative parameters for geyer cif.")
              if(par["sat"] > .Machine$integer.max-100)
		par["sat"] <- .Machine$integer.max-100
              par <- par[nms]
# Square r to avoid squaring in the Fortran code.
              par[[3]] <- par[[3]]**2
              return(par)
            },
            reach = function(par, ...) {
              r <- par[["r"]]
              g <- par[["gamma"]]
              return(if(g == 1) 0 else r)
            }
            ),
#       
# 9. The ``lookup'' device.  This permits simulating, at least
# approximately, ANY pairwise interaction function model
# with isotropic pair interaction (i.e. depending only on distance).
# The pair interaction function is provided as a vector of
# distances and corresponding function values which are used
# as a ``lookup table'' by the Fortran code.
#
       'lookup' = 
       list(
            fortran.id=9,
            multitype=FALSE,
            need.aux=FALSE,
            parhandler=function(par, ...) {
              nms <- c("beta","h")
              if(!is.list(par) || sum(!is.na(match(names(par),nms))) != 2) {
                stop(paste("For the lookup cif, par must be a named list\n",
                           "with components beta and h (and optionally r).\n",
                           "Bailing out."))
              }
              beta <- par[["beta"]]
              if(beta < 0)
		stop("Negative value of beta for lookup cif.")
              h.init <- par[["h"]]
              r <- par[["r"]]
              if(is.null(r)) {
		if(!is.stepfun(h.init))
                  stop(paste("For cif=lookup, if component r of",
                             "par is absent then component h must",
                             "be a stepfun object."))
		if(!is.cadlag(h.init))
                  stop(paste("The lookup pairwise interaction step",
			     "function must be right continuous,\n",
			     "i.e. built using the default values of the",
                             sQuote("f"), "and", sQuote("right"),
                             "arguments for stepfun."))
		r     <- knots(h.init)
		h0    <- get("yleft",envir=environment(h.init))
		h     <- h.init(r)
		nlook <- length(r)
		if(!identical(all.equal(h[nlook],1),TRUE))
                  stop(paste("The lookup interaction step function",
                             "must be equal to 1 for", dQuote("large"),
                             "distances."))
		if(r[1] <= 0)
                  stop(paste("The first jump point (knot) of the lookup",
                             "interaction step function must be",
                             "strictly positive."))
		h <- c(h0,h)
              } else {
		h     <- h.init
		nlook <- length(r)
		if(length(h) != nlook)
                  stop("Mismatch of lengths of h and r lookup vectors.")
		if(any(is.na(r)))
                  stop("Missing values not allowed in r lookup vector.")
		if(is.unsorted(r))
                  stop("The r lookup vector must be in increasing order.")
		if(r[1] <= 0)
                  stop(paste("The first entry of the lookup vector r",
                             "should be strictly positive."))
		h <- c(h,1)
              }
              if(any(h < 0))
		stop(paste("Negative values in the lookup",
                           "pairwise interaction function."))
              if(h[1] > 0 & any(h > 1))
		stop(paste("Lookup pairwise interaction function does",
                           "not define a valid point process."))
              rmax   <- r[nlook]
              r <- c(0,r)
              nlook <- nlook+1
              deltar <- mean(diff(r))
              if(identical(all.equal(diff(r),rep(deltar,nlook-1)),TRUE)) {
		equisp <- 1
		par <- c(beta,nlook,equisp,deltar,rmax,h)
              } else {
		equisp <- 0
		par <- c(beta,nlook,equisp,deltar,rmax,h,r)
               
              }
              return(par) 
            },
            reach = function(par, ...) {
              r <- par[["r"]]
              h <- par[["h"]]
              if(is.null(r)) 
                r <- knots(h)
              return(max(r))
            }
            )
  )
