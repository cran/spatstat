#
#    predict.ppm.S
#
#	$Revision: 1.48 $	$Date: 2008/08/12 07:24:04 $
#
#    predict.ppm()
#	   From fitted model obtained by ppm(),	
#	   evaluate the fitted trend or conditional intensity 
#	   at a grid/list of other locations 
#
#
# -------------------------------------------------------------------

predict.ppm <-
function(object, window, ngrid=NULL, locations=NULL,
         covariates=NULL, type="trend", 
         ..., check=TRUE, repair=TRUE) {
#
#  options for `type'
  type <- pickoption("type", type,
                     c(trend="trend", cif="cif", lambda="cif",
                       se="se", SE="se"))
#
#  guard against mistakes
#  
  dotargs <- list(...)
  nama <- names(dotargs)
  if(!is.null(nama) && any(nama == "newdata"))
    warning(paste("The use of the argument", sQuote("newdata"),
                  "is out-of-date. See help(predict.ppm)"))
  else if(length(dotargs) > 0)
    warning("Some arguments were ignored by predict.ppm")
#
#	'object' is the output of ppm()
#
  model <- object
  verifyclass(model, "ppm")
#  
  if(check && damaged.ppm(object)) {
    if(!repair)
      stop("object format corrupted; try update(object, use.internal=TRUE)")
    message("object format corrupted; repairing it.")
    object <- update(object, use.internal=TRUE)
  }
#
#       find out what kind of model it is
#
  mod <- summary(model, quick="no prediction")  # undocumented hack!
  stationary <- mod$stationary
  poisson    <- mod$poisson
  marked     <- mod$marked
  multitype  <- mod$multitype
  notrend    <- mod$no.trend
  trivial    <- poisson && notrend

  need.covariates <- mod$has.covars

  if(mod$antiquated)
    warning("The model was fitted by an out-of-date version of spatstat")
#
#       determine mark space
#  
  if(marked) {
    if(!multitype)
      stop("Prediction not yet implemented for general marked point processes")
    else 
      types <- levels(marks(mod$entries$data))
  }
#
#
#    Standard error only available for Poisson models
#
  if(type == "se" && !poisson)
      stop(paste("Standard error calculation",
                 "is only available for Poisson models"), call.=FALSE)
#
#      determine what kind of output is required:
#      (arguments present)    (output)  
#         window, ngrid    ->   image
#         locations (mask) ->   image
#         locations (other) ->  data frame
#
  if(!is.null(ngrid) && !is.null(locations))
    stop(paste("Only one of",
               sQuote("ngrid"), "and", sQuote("locations"),
               "should be specified"))

  if(is.null(ngrid) && is.null(locations)) 
    # use ngrid 
    ngrid <- 50
    
  want.image <- is.null(locations) ||
                    (is.owin(locations) && locations$type == "mask")
  make.grid <- !is.null(ngrid)
#
#
################   Determine prediction points  #####################
#  
  if(!want.image) {
    # (A) list of (x,y) coordinates given by `locations'
    xpredict <- locations$x
    ypredict <- locations$y
    if(is.null(xpredict) || is.null(ypredict)) {
      xy <- xy.coords(locations)
      xpredict <- xy$x
      xpredict <- xy$y
    }
    if(is.null(xpredict) || is.null(ypredict))
      stop(paste("Don't know how to extract x,y coordinates from",
                 sQuote("locations")))
    # marks if required
    if(marked)
      mpredict <- locations$marks  # locations is a data frame
  } else {
    # (B) pixel grid of points
    #
    if(!make.grid) 
    #    (B)(i) The grid is given in `locations'
      masque <- locations
    else {
    #    (B)(ii) We have to make the grid ourselves  
    #    Validate ngrid
    #  
      if(!is.null(ngrid)) {
        if(!is.numeric(ngrid))
          stop("ngrid should be a numeric vector")
        nn <- length(ngrid)
        if(nn < 1 || nn > 2)
          stop("ngrid should be a vector of length 1 or 2")
        if(nn == 1)
          ngrid <- rep(ngrid,2)
      }
      if(missing(window))
        window <- mod$entries$data$window
      masque <- as.mask(window, dimyx=ngrid)
    }
    # Hack -----------------------------------------------
    # gam with lo() will not allow extrapolation beyond the range of x,y
    # values actually used for the fit. Check this:
    tums <- termsinformula(model$trend)
    if(any(
           tums == "lo(x)" |
           tums == "lo(y)" |
           tums == "lo(x,y)" |
           tums == "lo(y,x)")
       ) {
      # determine range of x,y used for fit
      gg <- model$internal$glmdata
      gxr <- range(gg$x[gg$SUBSET])
      gyr <- range(gg$y[gg$SUBSET])
      # trim window to this range
      masque <- intersect.owin(masque, owin(gxr, gyr))
    }
    # ------------------------------------ End Hack
    #
    # Finally, determine x and y vectors for grid
    xx <- raster.x(masque)
    yy <- raster.y(masque)
    xpredict <- xx[masque$m]
    ypredict <- yy[masque$m]
  }

#
##################  CREATE DATA FRAME  ##########################
#                           ... to be passed to predict.glm()  
#
# First the x, y coordinates
  
    if(!marked) 
      newdata <- data.frame(x=xpredict, y=ypredict)
    else if(!want.image) 
      newdata <- data.frame(x=xpredict, y=ypredict, marks=mpredict)
    else {
      # replicate
      nt <- length(types)
      np <- length(xpredict)
      xpredict <- rep(xpredict,nt)
      ypredict <- rep(ypredict,nt)
      newdata <- data.frame(x = xpredict,
                            y = ypredict,
                            marks=rep(types, rep(np, nt)))
    }

#### Next the external covariates, if any
#
   if(need.covariates) {
     if(is.null(covariates)) {
       # Extract covariates from fitted model object
       # They have to be images.
       oldcov <- model$covariates
       if(is.null(oldcov))
         stop("External covariates are required, and are not available")
       if(is.data.frame(oldcov))
         stop(paste("External covariates are required.",
                    "Prediction is not possible at new locations"))
       covariates <- oldcov
     }
     covariates.df <-  mpl.get.covariates(covariates,
                       list(x=xpredict, y=ypredict), "prediction points")
     newdata <- cbind(newdata, covariates.df)
   }
#
######## Set up prediction variables ################################
#
#
# Provide SUBSET variable
#
        if(is.null(newdata$SUBSET))
          newdata$SUBSET <- rep(TRUE, nrow(newdata))
#
# Dig out information used in Berman-Turner device 
#        Vnames:     the names for the ``interaction variables''
#        glmdata:    the data frame used for the glm fit
#        glmfit:     the fitted glm object
#
  if(!trivial) {
	Vnames <- model$internal$Vnames
        glmdata <- getglmdata(model)
        glmfit <- getglmfit(model)
  }
  
############  COMPUTE PREDICTION ##############################
#
#   Compute the predicted value z[i] for each row of 'newdata'
#   Store in a vector z and reshape it later
#
###############################################################
  
  if(trivial) {
#############  UNIFORM POISSON PROCESS #####################

    lambda <- exp(model$coef[[1]])
    switch(type,
           cif=,
           trend={
             z <- rep(lambda, nrow(newdata))
           },
           se={
             npoints <- data.ppm(model)$n
             se.lambda <- lambda/sqrt(npoints)
             z <- rep(se.lambda, nrow(newdata))
           },
           stop("Internal error: unrecognised type"))
    
################################################################
  } else if((type %in% c("trend","se")) || poisson) {
#
#############  COMPUTE TREND ###################################
#	
#   set explanatory variables to zero
#	
    zeroes <- rep(0, nrow(newdata))    
    for(vn in Vnames)    
      newdata[[vn]] <- zeroes
#
#   invoke predict.glm()
#
    lambda <- predict(glmfit, newdata, type="response")
#
    switch(type,
           cif=,
           trend={
             z <- lambda
           },
           se={
             # extract variance-covariance matrix of parameters
             vc <- vcov(model)
             # compute model matrix
             mf <- model.frame(formula(model), newdata, ...)
             mm <- model.matrix(formula(model), mf, ...)
             # compute relative variance = diagonal of quadratic form
             nr <- nrow(newdata)
             vv <- numeric(nr)
             for(i in 1:nr) {
               mmi <- mm[i, ]
               vv[i] <- mmi %*% vc %*% mmi
             }
             z <- lambda * sqrt(vv)
           },
           stop("Internal error: unrecognised type"))

##############################################################  
  } else if(type == "cif" || type =="lambda") {
######### COMPUTE FITTED CONDITIONAL INTENSITY ################
#
# 	
  # set up arguments
    inter <- model$interaction
    X <- model$Q$data
    U <- list(x=newdata$x, y=newdata$y)
    Equal <- outer(X$x, U$x, "==") & outer(X$y, U$y, "==")
    if(marked) {
      U <- as.ppp(U, X$window)
      marks(U) <- Umarks <- newdata$marks
      Equal <- Equal & outer(marks(X), Umarks, "==")
    }
  # compute values of potential at the new sample points
    Vnew <- inter$family$eval(X, U, Equal,
                              inter$pot, inter$par, model$correction)
    if(!is.matrix(Vnew))
      stop("internal error: eval.pair.inter() did not return a matrix")

  # Negative infinite values signify cif = zero
    cif.equals.zero <- matrowany(Vnew == -Inf)
    
  # Insert the potential into the relevant column(s) of `newdata'
    if(ncol(Vnew) == 1)
      # Potential is real valued (Vnew is a column vector)
      # Assign values to a column of the same name in newdata
      newdata[[Vnames]] <- as.vector(Vnew)
      #
    else if(is.null(dimnames(Vnew)[[2]])) {
      # Potential is vector-valued (Vnew is a matrix)
      # with unnamed components.
      # Assign the components, in order of their appearance,
      # to the columns of newdata labelled Vnames[1], Vnames[2],... 
      for(i in seq(Vnames))
        newdata[[Vnames[i] ]] <- Vnew[,i]
      #
    } else {
      # Potential is vector-valued (Vnew is a matrix)
      # with named components.
      # Match variables by name
      for(vn in Vnames)    
        newdata[[vn]] <- Vnew[,vn]
      #
    }
  # invoke predict.glm
  z <- predict(glmfit, newdata, type="response")

  # reset to zero if potential was zero
  if(any(cif.equals.zero))
    z[cif.equals.zero] <- 0
    
#################################################################    
  } else
     stop(paste("Unrecognised type", sQuote(type)))

#################################################################
#
# reshape the result
#
    if(!want.image) 
      out <- as.vector(z)
    else {
      # make an image of the right shape
      imago <- as.im(masque)
      imago <- eval.im(as.double(imago))
      if(!marked) {
        # single image
        out <- imago
        # set entries
        out$v[masque$m] <- z
      } else {
        # list of images
        out <- list()
        for(i in seq(types)) {
          outi <- imago
          # set entries
          outi$v[masque$m] <- z[newdata$marks == types[i]]
          out[[i]] <- outi
        }
        names(out) <- paste("mark",types, sep="")
        class(out) <- c("listof", class(out))
      }
    }
#  
#  FINISHED
#  
  return(out)
}


####################################################################
#
# compute pointwise uncertainty of fitted intensity
#
model.se.image <- function(fit, W=as.owin(fit), ..., what="sd") {
  if(!is.poisson.ppm(fit))
    stop("Only implemented for Poisson point process models", call.=FALSE)
  what <- pickoption("option", what,
                     c(sd="sd", var="var", cv="cv", CV="cv", ce="ce", CE="ce"))
  W <- as.mask(as.owin(W))
  # variance-covariance matrix of coefficients
  vc <- vcov(fit)
  np <- dim(vc)[1]
  # extract sufficient statistic for each coefficient
  mm <- model.images(fit, W, ...)
  # compute fitted intensity 
  lam <- predict(fit, locations=W)
  # initialise resulting image
  U <- as.im(W)
  U[] <- 0
  # compute pointwise matrix product, assuming vc is symmetric
  for(i in 1:np) {
    Si <- mm[[i]]
    aii <- vc[i,i]
    U <- eval.im(U + aii * Si^2)
    if(i > 1) {
      for(j in 1:(i-1)) {
        Sj <- mm[[j]]
        aij <- vc[i,j]
        twoaij <- 2 * aij
        U <- eval.im(U + twoaij * Si * Sj)
      }
    }
  }
  # the matrix product is the relative variance (CV)
  if(what=="cv")
    return(U)
  # relative sd
  if(what=="ce") {
    U <- eval.im(sqrt(U))
    return(U)
  }
  # multiply by squared intensity to obtain variance
  U <- eval.im(U * lam^2)
  # variance
  if(what=="var")
    return(U)
  # compute SD and return
  U <- eval.im(sqrt(U))
  return(U)
}
