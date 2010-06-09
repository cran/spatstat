#
#  rmhmodel.ppm.R
#
#   convert ppm object into format palatable to rmh.default
#
#  $Revision: 2.37 $   $Date: 2010/06/05 03:41:28 $
#
#   .Spatstat.rmhinfo
#   rmhmodel.ppm()
#

.Spatstat.Rmhinfo <-
list(
     "Lennard-Jones process" =
     function(coeffs, inte) {
       sigma   <- inte$par$sigma
       epsilon <- inte$par$epsilon
       return(list(cif='lennard',
                   par=list(sigma=sigma, epsilon=epsilon)))
     },
     "Fiksel process" =
     function(coeffs, inte) {
       hc <- inte$par$hc
       r  <- inte$par$r
       kappa <- inte$par$kappa
       a <- inte$interpret(coeffs,inte)$param$a
       return(list(cif='fiksel',
                   par=list(r=r,hc=hc,kappa=kappa,a=a)))
     },
     "Diggle-Gratton process" =
     function(coeffs, inte) {
       kappa <- inte$interpret(coeffs,inte)$param$kappa
       delta <- inte$par$delta
       rho   <- inte$par$rho
       return(list(cif='diggra',
                   par=list(kappa=kappa,delta=delta,rho=rho)))
     },
     "Hard core process" =
     function(coeffs, inte) {
       hc <- inte$par$hc
       return(list(cif='hardcore',
                   par=list(hc=hc)))
     },
     "Geyer saturation process" =
     function(coeffs, inte) {
       gamma <- inte$interpret(coeffs,inte)$param$gamma
       r <- inte$par$r
       sat <- inte$par$sat
       return(list(cif='geyer',
                   par=list(gamma=gamma,r=r,sat=sat)))
     },
     "Soft core process" =
     function(coeffs, inte) {
       kappa <- inte$par$kappa
       sigma <- inte$interpret(coeffs,inte)$param$sigma
       return(list(cif="sftcr",
                   par=list(sigma=sigma,kappa=kappa)))
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
                   par=list(gamma=gamma,r=r,hc=hc)))
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
                   par=list(gamma=gamma,iradii=iradii,hradii=hradii)))
     },
     "Piecewise constant pairwise interaction process" =
     function(coeffs, inte) {
       r <- inte$par$r
       gamma <- (inte$interpret)(coeffs, inte)$param$gammas
       h <- stepfun(r, c(gamma, 1))
       return(list(cif='lookup', par=list(h=h)))
     },
     "Area-interaction process" =
     function(coeffs, inte) {
       r <- inte$par$r
       eta <- (inte$interpret)(coeffs, inte)$param$eta
       return(list(cif='areaint', par=list(eta=eta,r=r)))
     },
     "hybrid Geyer process" =
     function(coeffs, inte) {
       r <- inte$par$r
       sat <- inte$par$sat
       gamma <- (inte$interpret)(coeffs,inte)$param$gammas
       return(list(cif='badgey',par=list(gamma=gamma,r=r,sat=sat)))
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
#           OrdThresh                <none>
#
#  'dgs' has no canonical parameters (it's determined by its irregular
#  parameter rho) so there can't be an interaction object for it.
#
#  Implementing rmh.default for the others is probably too hard.


rmhmodel.ppm <- function(model, win, ..., verbose=TRUE, project=TRUE,
                         control=rmhcontrol()) {
  # converts ppm object `model' into format palatable to rmh.default
  
    verifyclass(model, "ppm")
    X <- model

    if(verbose)
      cat("Extracting model information...")
    
    # Extract essential information
    Y <- summary(X)

    if(Y$marked && !Y$multitype)
      stop("Not implemented for marked point processes other than multitype")

    if(Y$uses.covars && is.data.frame(X$covariates))
      stop(paste("This model cannot be simulated, because the",
                 "covariate values were given as a data frame."))
    
    # enforce defaults for `control'

    control <- rmhcontrol(control)
    
    ########  Interpoint interaction
    if(Y$poisson) {
      Z <- list(cif="poisson",
                par=list())  # par is filled in later
    } else {
      # First check version number of ppm object
      if(Y$antiquated) 
        stop(paste("This model was fitted by a very old version",
                   "of the package: spatstat", Y$version,
                   "; simulation is not possible.",
                   "Re-fit the model using your original code"))
      else if(Y$old)
        warning(paste("This model was fitted by an old version",
                      "of the package: spatstat", Y$version,
                      ". Re-fit the model using update.ppm",
                      "or your original code"))
      # Extract the interpoint interaction object
      inte <- Y$entries$interaction
      # Determine whether the model can be simulated using rmh
      siminfo <- .Spatstat.Rmhinfo[[inte$name]]
      if(is.null(siminfo))
        stop(paste("Simulation of a fitted", sQuote(inte$name),
                   "has not yet been implemented"))
      
      # Get fitted model's canonical coefficients
      coeffs <- Y$entries$coef
      if(newstyle.coeff.handling(inte)) {
        # extract only the interaction coefficients
        Vnames <- Y$entries$Vnames
        coeffs <- coeffs[Vnames]
      }
      # Ensure the fitted model is valid
      # (i.e. exists mathematically as a point process)
      if(!valid.ppm(model)) {
        if(project) {
          if(verbose)
            cat("Model is invalid - projecting it\n")
          if(is.null(inte$project))
            stop("Internal error: interaction has no projection operator")
          coeffs <- inte$project(coeffs, inte)
        }
        else
          stop("The fitted model is not a valid point process")
      }
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
       
    ######## Window for result 
    
    if(missing(win))
      win <- Y$entries$data$window

    Z$w <- win

    ######## Expanded window for simulation?

    covims <- if(Y$uses.covars) X$covariates[Y$covars.used] else NULL
    
    wsim <- rmhResolveExpansion(win, control, covims, "covariate")$wsim
      
    ###### Trend or Intensity ############

    if(verbose)
      cat("Evaluating trend...")
    
    if(Y$stationary) {
      # first order terms (beta or beta[i]) are carried in Z$par$beta
      Z$par[["beta"]] <- Y$trend$value
      # Note the construction m[["name"]] <- value
      # ensures that m does not need to be coerced from vector to list
      Z$trend <- NULL
    } else {
      # trend terms present
      # all first order effects are subsumed in Z$trend
      Z$par[["beta"]] <- if(!Y$marked) 1 else rep(1, length(Z$types))
      # predict on window possibly larger than original data window
      Z$trend <- 
        if(control$condtype != "none" && wsim$type == "mask")
          predict(X, window=wsim, type="trend", locations=wsim)
        else 
          predict(X, window=wsim, type="trend")
    }
    if(verbose)
      cat("done.\n")
    Z <- rmhmodel(Z)
    return(Z)
}

rmhResolveExpansion <- function(win, control, imagelist, itype="covariate") {
  # Determine expansion window for simulation

  if(control$force.noexp) {
    # Expansion prohibited
    return(list(wsim=win, expanded=FALSE))
  }
  
  # Determine proposed expansion window
  expand <- control$expand
  if(is.null(expand)) {
    want.expand <- FALSE
    wexp <- win
  } else if(is.owin(expand)) {
    want.expand <- TRUE
    wexp <- expand
  } else if(is.numeric(expand) && length(expand) == 1) {
    if(expand < 1)
      stop("expand must be >= 1")
    wexp <- if(expand == 1) win else expand.owin(win, expand)
    want.expand <- (expand > 1)
  } else stop(paste("Unrecognised format for", sQuote("expand")))
    
  # Decide whether to expand

  if(!want.expand)
    return(list(wsim=win, expanded=FALSE))

  isim <- unlist(lapply(imagelist, is.im))
  imagelist <- imagelist[isim]

  if(length(imagelist) == 0) {
    # Unlimited expansion is feasible
    return(list(wsim=wexp, expanded=TRUE))
  }

  # Expansion is limited to domain of image data
  # Determine maximum possible expansion window
  wins <- lapply(imagelist, as.owin)
  if(length(wins) == 1) {
    cwin <- wins[[1]]
  } else {
    names(wins) <- NULL
    cwin <- do.call("intersect.owin", wins)
  }
  
  if(!is.subset.owin(wexp, cwin)) {
    # Cannot expand to proposed window
    if(control$force.exp)
      stop(paste("Cannot expand the simulation window,",
                 "because the", itype, "images do not cover",
                 "the expanded window"))
      # Take largest possible window
    wexp <- intersect.owin(wexp, cwin)
  }
  return(list(wsim=wexp, expanded=TRUE))
}

