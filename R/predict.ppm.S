#
#    predict.ppm.S
#
#	$Revision: 1.15 $	$Date: 2004/01/08 13:12:27 $
#
#    predict.ppm()
#	   From fitted model obtained by mpl(),	
#	   evaluate the fitted trend or conditional intensity 
#	   at a grid/list of other locations 
#
#
# -------------------------------------------------------------------

predict.ppm <-
function(object, newdata, window, ngrid=NULL, type="trend", ...) {
#
#	'object' is the output of mpl()
#
  model <- object
  verifyclass(model, "ppm")
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
#
#       determine mark space
#  
  if(marked) {
    if(!multitype)
      stop("Not yet implemented for general marked point processes")
    else 
      types <- levels(mod$entries$data$marks)
  }
#
#      determine what kind of output is required: image(s) or data frame
#  
  want.image <- missing(newdata)
#
#  Create data frame to pass to predict.glm()
#
  if(want.image) {
#
#      validate parameters for image
#  
    if(missing(window))
      window <- mod$entries$data$window
    if(!is.null(ngrid)) {
      if(!is.numeric(ngrid))
        stop("ngrid should be a numeric vector")
      nn <- length(ngrid)
      if(nn < 1 || nn > 2)
        stop("ngrid should be a vector of length 1 or 2")
      if(nn == 1)
        ngrid <- rep(ngrid,2)
    }
    masque <- as.mask(window, dimyx=ngrid)
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
  # determine x and y vectors for grid
    xx <- raster.x(masque)
    yy <- raster.y(masque)
    xx <- xx[masque$m]
    yy <- yy[masque$m]
    if(!marked) 
      newdata <- data.frame(x=xx, y=yy)
    else {
      nt <- length(types)
      nx <- length(xx)
      newdata <- data.frame(x = rep(xx, nt),
                            y = rep(yy, nt),
                            marks=rep(types, rep(nx, nt)))
    }
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
# Dig out information used in the original call to mpl().
#        Vnames:     the names for the ``interaction variables''
#        glmdata:    the data frame used for the glm fit
# and set the contrasts according to the convention followed in mpl().
#
  if(!trivial) {
	Vnames <- model$internal$Vnames

        if(exists("is.R") && is.R())
	    glmdata <- model$internal$glmdata
        else        
            assign("glmdata", model$internal$glmdata, f = 1)
  }
  
############  COMPUTE PREDICTION ##############################
#
#   Compute the predicted value z[i] for each row of 'newdata'
#   Store in a vector z and reshape it later
#

###############################################################  
  if(trivial) {
#############  COMPUTE CONSTANT INTENSITY #####################

    lambda <- exp(model$theta[[1]])
    z <- rep(lambda, nrow(newdata))
    
################################################################
  } else if(type == "trend" || poisson) {
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
    z <- predict(model$internal$glmfit, newdata, type="response")

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
      U$marks <- newdata$marks
      Equal <- Equal & outer(X$marks, U$marks, "==")
    }
  # compute values of potential at the new sample points
    Vnew <- inter$family$eval(X, U, Equal,
                              inter$pot, inter$par, model$correction)
    if(!is.matrix(Vnew))
      stop("internal error: eval.pair.inter() did not return a matrix")
  
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
  z <- predict(model$internal$glmfit, newdata, type="response")

#################################################################    
  } else
     stop(paste("Unrecognised type \'", type, "\'\n", sep=""))

#################################################################
#
# reshape the result
#
    if(!want.image) 
      out <- as.vector(z)
    else {
      # make an image of the right shape
      imago <- as.im(masque)
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
      }
    }

####################################################################
#
#  
  invisible(out)
}

