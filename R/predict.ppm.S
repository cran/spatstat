#
#    predictppm.S
#
#	$Revision: 1.11 $	$Date: 2002/01/18 06:57:40 $
#
#    predict.ppm()
#	   From fitted model obtained by mpl(),	
#	   evaluate the fitted trend or conditional intensity 
#	   at a grid/list of other locations 
#
#
# -------------------------------------------------------------------

predict.ppm <-
function(object, newdata, nx = 40, ny = NULL, type="trend", ...) {
#
#	'object' is the output of mpl()
#
  model <- object
  verifyclass(model, "ppm")
#
#       find out what kind of model it is
#
  stationary <- is.stationary.ppm(model)
  poisson <- is.poisson.ppm(model)
  markeddata <- is.marked(model$Q$data)
  notrend <- no.trend.ppm(model)
  trivial <- poisson && notrend

  if(markeddata) {
      mrks <- model$Q$data$marks
  }    
#
#  Create data frame to pass to predict.glm()
#		
  want.grid <- missing(newdata)

  if(want.grid) {
  # create a rectangular grid of nx x ny points in window
    window <- model$Q$data$window
    if(is.null(ny))
      ny <- nx
    xr <- window$xrange
    yr <- window$yrange
  # gam with lo() will not allow extrapolation beyond the range of x,y
  # values actually used for the fit. Check this:
    tums <- termsinformula(model$trend)
    if(any(
           tums == "lo(x)" |
           tums == "lo(y)" |
           tums == "lo(x,y)" |
           tums == "lo(y,x)")
    ) {
      gg <- model$internal$glmdata
      gxr <- range(gg$x[gg$SUBSET])
      gyr <- range(gg$y[gg$SUBSET])
      intersect.ranges <- function(r, s) {
        c(max(r[1],s[1]), min(r[2],s[2]))
      }
      xr <- intersect.ranges(xr, gxr)
      yr <- intersect.ranges(yr, gyr)
    }
  # determine x and y vectors for nx * ny grid
    dx <- diff(xr)/nx
    dy <- diff(yr)/ny
    xvals <- seq(xr[1] + dx/2, xr[2]-dx/2, length=nx)
    yvals <- seq(yr[1] + dy/2, yr[2]-dy/2, length=ny)
    if(!markeddata) {
      # create nx * ny grid
      newdata <- expand.grid(x=xvals, y = yvals)
      # initialise output list
      out <- list(x=xvals, y=yvals)
    } else {
      if(!is.factor(mrks))
        stop("marks are not a factor")
      mrkvals <- levels(mrks)
      # create (nx * ny * nmarks) marked points
      newdata <- expand.grid(x=xvals, y=yvals, marks=mrkvals)
      # initialise output list
      out <- list(x=xvals, y=yvals, marks = mrkvals)
    }
  }

#
######## Set up prediction variables ################################
#
#
# Provide SUBSET variable
#
        if(is.null(newdata$SUBSET))
          newdata$SUBSET <- rep(T, nrow(newdata))
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
    if(markeddata) {
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
    if(!want.grid) 
      out <- as.vector(z)
    else {
      out$z <-
        if(!markeddata)
          matrix(z, nx, ny)
        else
          array(z, dim=c(nx,ny,length(mrkvals)),
                dimnames=list(NULL, NULL, as.character(mrkvals)))
    }

####################################################################
#
#  
  invisible(out)
}

