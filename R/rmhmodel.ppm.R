#
#  rmhmodel.ppm.R
#
#   convert ppm object into format palatable to rmh.default
#
#  $Revision: 2.4 $   $Date: 2004/01/27 11:28:21 $
#
#   .Spatstat.rmhinfo
#   rmhmodel.ppm()
#

.Spatstat.Rmhinfo <-
list(
     "Diggle-Gratton process" =
     function(coeffs, inte) {
       kappa <- inte$interpret(coeffs,inte)$param$kappa
       delta <- inte$par$delta
       rho   <- inte$par$rho
       return(list(cif='diggra',
                   par=c(kappa=kappa,delta=delta,rho=rho)))
     },
     "Geyer saturation process" =
     function(coeffs, inte) {
       gamma <- inte$interpret(coeffs,inte)$param$gamma
       r <- inte$par$r
       sat <- inte$par$saturate
       return(list(cif='geyer',
                   par=c(gamma=gamma,r=r,sat=sat)))
     },
     "Soft core process" =
     function(coeffs, inte) {
       kappa <- inte$par$kappa
       sigma <- inte$interpret(coeffs,inte)$param$sigma
       return(list(cif="sftcr",
                   par=c(sigma=sigma,kappa=kappa)))
     },
     "Strauss process" =
     function(coeffs, inte) {
       gamma <- inte$interpret(coeffs,inte)$param$gamma
       r <- inte$par$r
       return(list(cif = "strauss",
                   par = list(gamma = gamma, r = r)))
     },
     "Strauss - hard core process" =
     function(coeffs, inte) {
       gamma <- inte$interpret(coeffs,inte)$param$gamma
       r <- inte$par$r
       hc <- inte$par$hc
       return(list(cif='straush',
                   par=c(gamma=gamma,r=r,hc=hc)))
     },
     "Multitype Strauss process" =
     function(coeffs, inte) {
       # interaction radii r[i,j]
       radii <- inte$par$radii
       # interaction parameters gamma[i,j]
       gamma <- (inte$interpret)(coeffs, inte)$param$gammas
       return(list(cif='straussm',
                   par=list(gamma=gamma,radii=radii)))
     },
     "Multitype Strauss Hardcore process" =
     function(coeffs, inte) {
       # interaction radii r[i,j]
       iradii <- inte$par$iradii
       # hard core radii r[i,j]
       hradii <- inte$par$hradii
       # interaction parameters gamma[i,j]
       gamma <- (inte$interpret)(coeffs, inte)$param$gammas
       return(list(cif='straushm',
                   par=list(gamma=gamma,iradii=radii,hradii=hradii)))
     }
)


# OTHER MODELS not yet implemented:
#
#
#      interaction object           rmh.default 
#      ------------------           -----------
#
#           <none>                   dgs
#
#           PairPiece                <none>
#
#           LennardJones             <none>
#
#           OrdThresh                <none>
#
#  'dgs' has no canonical parameters (it's determined by its irregular
#  parameter rho) so there can't be an interaction object for it.
#
#  Implementing rmh.default for the others is probably too hard.


rmhmodel.ppm <- function(model, win, ..., verbose=TRUE, project=TRUE) {
  # convert ppm object `model' into format palatable to rmh.default
  
    verifyclass(model, "ppm")
    X <- model

    # Extract essential information
    Y <- summary(X)

    if(Y$marked && !Y$multitype)
      stop("Not implemented for marked point processes other than multitype")

    ########  Interpoint interaction
    if(Y$poisson) {
      Z <- list(cif="poisson",
                par=list())  # par is filled in later
    } else {
      # First check version number of ppm object
      ver <- model$version
      if(is.null(ver)
         || is.character(ver)   # old style, version <= 1.3-4
         || is.null(major <- ver$major)
         || is.null(minor <- ver$minor)
         || (major == 1 && minor < 4)) {
        whinge <- paste(
        "This model was fitted by an earlier version of spatstat;\n",
        "simulation is not possible.\n",
        "Re-fit the model using the current version of the package")
        stop(whinge)
      }
      # Extract the interpoint interaction object
      inte <- Y$entries$interaction
      # Get fitted model's canonical coefficients
      coeffs <- Y$entries$theta
      # Ensure the fitted model is valid
      # (i.e. exists mathematically as a point process)
      valid <- inte$valid(coeffs, inte)
      # if not, 
      if(!valid) {
        if(project) {
          if(verbose)
            cat("Model is invalid - projecting it\n")
          coeffs <- inte$project(coeffs, inte)
        }
        else
          stop("The fitted model is not a valid point process")
      }
      # Determine whether the model can be simulated using rmh
      siminfo <- .Spatstat.Rmhinfo[[inte$name]]
      if(is.null(siminfo))
        stop(paste("Simulation of a fitted \'", inte$name,
                   "\' has not yet been implemented"))
      # Translate the model to the format required by rmh.default
      Z <- siminfo(coeffs, inte)
      if(is.null(Z))
        stop("The model cannot be simulated")
      else if(is.null(Z$cif))
        stop(paste("Internal error: no cif returned from .Spatstat.Rmhinfo"))
    }

    # Don't forget the types
    if(Y$multitype && is.null(Z$types))
      Z$types <- levels(Y$entries$marks)
       
    ######## Simulation window 
    
    if(missing(win))
      win <- Y$entries$data$window

    Z$w <- win

    ###### Trend or Intensity ############

    if(Y$stationary) {
      # first order terms (beta or beta[i]) are carried in Z$par$beta
      Z$par[["beta"]] <- Y$trend$value
      Z$trend <- NULL
    } else {
      # all first order effects are subsumed in Z$trend
      Z$par[["beta"]] <- if(!Y$marked) 1 else rep(1, length(Z$types))
      Z$trend <- predict(X, window=as.mask(win), type="trend")
    }
    # Note the construction m[["name"]] <-
    # ensures that m does not need to be coerced from vector to list
    return(Z)
}

