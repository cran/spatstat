#
#	$Revision: 4.10 $	$Date: 2003/03/11 03:24:03 $
#
#    mpl()
#          Fit a point process model to a two-dimensional point pattern
#          The model may include
#                  - trend (arbitrary);
#                  - dependence on covariates;  
#                   - arbitrary interaction structure
#			(with p-dimensional regular parameter)
#	   The window of observation may have arbitrary shape.
#
#
# -------------------------------------------------------------------
#	
"mpl" <- 
function(Q,
         trend = ~1,
	 interaction = NULL,
         data,
	 correction="border",
	 rbord = 0,
         use.gam=FALSE
) {
#
# Extract quadrature scheme 
#
	if(verifyclass(Q, "ppp", fatal = FALSE)) {
#		warning("using default quadrature scheme")
		Q <- quadscheme(Q)   
	} else if(!verifyclass(Q, "quad", fatal=FALSE))
		stop("First argument Q should be a quadrature scheme")
#
# Data points
  X <- Q$data
#
# Data and dummy points together 
  P <- union.quad(Q)
#
#
# Interpret the call
want.trend <- !is.null(trend) && !identical.formulae(trend, ~1)
want.inter <- !is.null(interaction) && !is.null(interaction$family)

the.call <- deparse(sys.call())
the.version <- "Spatstat 1.3; mpl.S $Date: 2003/03/11 03:24:03 $"

if(use.gam && exists("is.R") && is.R()) 
  require(mgcv)
        
if(!want.trend && !want.inter) {
  # the model is the uniform Poisson process
  # The MPLE of its intensity is the MLE
  npts <- X$n
  volume <- area.owin(X$window) * markspace.integral(X)
  lambda <- npts/volume
  theta <- list("log(lambda)"=log(lambda))
  maxlogpl <- npts * (log(lambda) - 1)
  rslt <- list(
               theta       = theta,
               coef        = theta,
               trend       = NULL,
               interaction = NULL,
               Q           = Q,
               maxlogpl    = maxlogpl,
               internal    = list(),
	       correction  = correction,
               rbord       = rbord,
               call        = the.call,
               version     = the.version)
  class(rslt) <- "ppm"
  return(rslt)
}


################################################################
################ C o m p u t e     d a t a  ####################
################################################################


### Form the weights and the ``response variable''.

W <- w.quad(Q)
Z <- is.data(Q)
Y <- Z/W
n <- length(Y)

MARKS <- marks.quad(Q)    # is NULL for unmarked patterns
	
glmdata <- data.frame(W=W, Y=Y)
SUBSET <- rep(TRUE, n)

zeroes <- attr(W, "zeroes")
if(!is.null(zeroes))
	SUBSET <-  !zeroes

####################### T r e n d ##############################

if(want.trend) {
  # Default explanatory variables for trend
  glmdata <- data.frame(glmdata, x=P$x, y=P$y)
  if(!is.null(MARKS))
    glmdata <- data.frame(glmdata, marks=MARKS)
  # 
  if(!missing(data)) {
# Append `data' to `glmdata'
# First validate 'data'
#   Check it has the right number of rows
    if(nrow(data) != nrow(glmdata))
      stop("Number of rows of \"data\" != number of points in \"Q\"")
#   Check any duplication of reserved names
    name.match <- outer(names(glmdata), names(data), "==")
    if(any(name.match)) {
      is.matched <- apply(name.match, 2, any)
      matched.names <- (names(data))[is.matched]
      if(sum(is.matched) == 1) {
        stop(paste("the variable name \"",
                   matched.names,
                   "\" in \'data\' is a reserved name - see help(mpl)"))
      } else {
        stop(paste("the variable names \"",
                   paste(matched.names, collapse=", "),
                   "\" in \'data\' are reserved names - see help(mpl)"))
      }
    }
#   OK, append 'data'    
    glmdata <- data.frame(glmdata,data)
  }
}

###################### I n t e r a c t i o n ####################

Vnames <- NULL

if(want.inter) {

  verifyclass(interaction, "interact")
  
  # Calculations require a matrix (data) x (data + dummy) indicating equality
  E <- equals.quad(Q)
  
  # Form the matrix of "regression variables" V.
  # The rows of V correspond to the rows of P (quadrature points)
  # while the column(s) of V are the regression variables (log-potentials)

  V <- interaction$family$eval(X, P, E,
                        interaction$pot,
                        interaction$par,
                        correction)

  if(!is.matrix(V))
    stop("interaction evaluator did not return a matrix")

  # Augment data frame by appending the regression variables for interactions.
  #
  # If there are no names provided for the columns of V,
  # call them "Interact.1", "Interact.2", ...

  if(is.null(dimnames(V)[[2]])) {
    # default names
    nc <- ncol(V)
    dimnames(V) <- list(dimnames(V)[[1]], 
      if(nc == 1) "Interaction" else paste("Interact.", 1:nc, sep=""))
  }
  Vnames <- dimnames(V)[[2]]
  
  glmdata <- data.frame(glmdata, V)   

# Keep only those quadrature points for which the
# conditional intensity is nonzero. 

#KEEP  <- apply(V != -Inf, 1, all)
KEEP  <- matrowall(V != -Inf)

SUBSET <- SUBSET & KEEP

if(any(Z & !KEEP)) {
        howmany <- sum(Z & !KEEP)
	warning(paste(howmany, "data point(s) are illegal (zero conditional intensity under the model)"))
#	browser()
}

}
######################################################################
################ F I T    M O D E L  #################################
######################################################################

# Determine the domain of integration for the pseudolikelihood.

if(correction == "border" && !missing(rbord)) {
	bd <- bdist.points(P)
	DOMAIN <- (bd >= rbord)
	SUBSET <- DOMAIN & SUBSET
}

glmdata <- data.frame(glmdata, SUBSET=SUBSET)

#################  F o r m u l a   ##################################


if(!want.trend) trend <- ~1 
trendpart <- paste(as.character(trend), collapse=" ")
rhs <- paste(c(trendpart, Vnames), collapse= "+")
fmla <- paste("Y ", rhs)
fmla <- as.formula(fmla)

################# F i t    i t   ####################################

# Fit the generalized linear/additive model.

if(want.trend && use.gam)
  FIT  <- gam(fmla, family=quasi(link=log, var=mu), weights=W,
              data=glmdata, subset=(SUBSET=="TRUE"), maxit=50)
else
  FIT  <- glm(fmla, family=quasi(link=log, var=mu), weights=W,
              data=glmdata, subset=(SUBSET=="TRUE"), maxit=50)
  
################  I n t e r p r e t    f i t   #######################

# Fitted coefficients

co <- FIT$coef
theta <- if(exists("is.R") && is.R()) NULL else dummy.coef(FIT)

  
# attained value of max log pseudolikelihood
maxlogpl <-  -(deviance(FIT)/2 + sum(log(W[Z & SUBSET])) + sum(Z & SUBSET))

######################################################################
# Clean up & return 

rslt <- list(
             theta        = theta,
             coef         = co,
             trend        = if(want.trend) trend       else NULL,
             interaction  = if(want.inter) interaction else NULL,
             Q            = Q,
             maxlogpl     = maxlogpl, 
             internal     = list(glmfit=FIT, glmdata=glmdata, Vnames=Vnames),
             correction   = correction,
             rbord        = rbord,
             call         = the.call,
             version      = the.version)
class(rslt) <- "ppm"
return(rslt)
}  

