#    mpl.R
#
#	$Revision: 5.13 $	$Date: 2004/09/01 04:13:03 $
#
#    mpl.engine()
#          Fit a point process model to a two-dimensional point pattern
#          by maximum pseudolikelihood
#
#    mpl.prepare()
#          set up data for glm procedure
#
# -------------------------------------------------------------------
#

"mpl" <- function(Q,
         trend = ~1,
	 interaction = NULL,
         data = NULL,
	 correction="border",
	 rbord = 0,
         use.gam=FALSE) {
   .Deprecated("ppm", package="spatstat")
   ppm(Q, trend, interaction, data, correction, rbord, use.gam, method="mpl")
}

"mpl.engine" <- 
function(Q,
         trend = ~1,
	 interaction = NULL,
         covariates = NULL,
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

the.version <- list(major=1,
                    minor=5,
                    release=3,
                    date="$Date: 2004/09/01 04:13:03 $")

if(use.gam && exists("is.R") && is.R()) 
  require(mgcv)
        
if(!want.trend && !want.inter) {
  # the model is the uniform Poisson process
  # The MPLE (= MLE) can be evaluated directly
  npts <- X$n
  volume <- area.owin(X$window) * markspace.integral(X)
  lambda <- npts/volume
  theta <- list("log(lambda)"=log(lambda))
  maxlogpl <- npts * (log(lambda) - 1)
  rslt <- list(
               method      = "mpl",
               theta       = theta,
               coef        = theta,
               trend       = NULL,
               interaction = NULL,
               Q           = Q,
               maxlogpl    = maxlogpl,
               internal    = list(),
	       correction  = correction,
               rbord       = rbord,
               version     = the.version)
  class(rslt) <- "ppm"
  return(rslt)
}

        
#################  P r e p a r e    D a t a   ######################
        
prep <- mpl.prepare(Q, X, P, trend, interaction,
                    covariates, want.trend, want.inter, correction, rbord)

fmla <- prep$fmla
glmdata <- prep$glmdata

        
################# F i t    i t   ####################################

# Fit the generalized linear/additive model.

if(want.trend && use.gam)
  FIT  <- gam(fmla, family=quasi(link=log, var=mu), weights=.mpl.W,
              data=glmdata, subset=(.mpl.SUBSET=="TRUE"),
              control=gam.control(maxit=50))
else
  FIT  <- glm(fmla, family=quasi(link=log, var=mu), weights=.mpl.W,
              data=glmdata, subset=(.mpl.SUBSET=="TRUE"),
              control=glm.control(maxit=50))
  
################  I n t e r p r e t    f i t   #######################

# Fitted coefficients

co <- FIT$coef
theta <- if(exists("is.R") && is.R()) NULL else dummy.coef(FIT)

     W <- glmdata$.mpl.W
SUBSET <- glmdata$.mpl.SUBSET        
     Z <- is.data(Q)
Vnames <- prep$Vnames
        
# attained value of max log pseudolikelihood
maxlogpl <-  -(deviance(FIT)/2 + sum(log(W[Z & SUBSET])) + sum(Z & SUBSET))

######################################################################
# Clean up & return 

rslt <- list(
             method       = "mpl",
             theta        = theta,
             coef         = co,
             trend        = if(want.trend) trend       else NULL,
             interaction  = if(want.inter) interaction else NULL,
             Q            = Q,
             maxlogpl     = maxlogpl, 
             internal     = list(glmfit=FIT, glmdata=glmdata, Vnames=Vnames),
             covariates   = covariates,
             correction   = correction,
             rbord        = rbord,
             version      = the.version)
class(rslt) <- "ppm"
return(rslt)
}  


##########################################################################
### /////////////////////////////////////////////////////////////////////
##########################################################################


mpl.prepare <- function(Q, X, P, trend, interaction, covariates, 
                        want.trend, want.inter, correction, rbord) {

# Validate/evaluate covariates
if(want.trend && !is.null(covariates))
  covariates.df <- mpl.get.covariates(covariates, P, "quadrature points")

################ C o m p u t e     d a t a  ####################

        
### Form the weights and the ``response variable''.

.mpl <- list()
.mpl$W <- w.quad(Q)
.mpl$Z <- is.data(Q)
.mpl$Y <- .mpl$Z/.mpl$W
.mpl$MARKS <- marks.quad(Q)  # is NULL for unmarked patterns

glmdata <- data.frame(.mpl.W = .mpl$W,
                      .mpl.Y = .mpl$Y)
        
n <- nrow(glmdata)
.mpl$SUBSET <- rep(TRUE, n)
	
internal.names <- c(".mpl.W", ".mpl.Y", ".mpl.Z", ".mpl.SUBSET",
                    "SUBSET", ".mpl")

reserved.names <- c("x", "y", "marks", internal.names)
                    
zeroes <- attr(.mpl$W, "zeroes")
if(!is.null(zeroes))
	.mpl$SUBSET <-  !zeroes

####################### T r e n d ##############################

  check.clashes <- function(forbidden, offered, where) {
    name.match <- outer(forbidden, offered, "==")
    if(any(name.match)) {
      is.matched <- apply(name.match, 2, any)
      matched.names <- (offered)[is.matched]
      if(sum(is.matched) == 1) {
        return(paste("The variable \"",
                   matched.names,
                   "\" in ", where,
                   " is a reserved name", sep=""))
      } else {
        return(paste("The variables \"",
                   paste(matched.names, collapse="\", \""),
                   "\" in ", where,
                   " are reserved names", sep=""))
      }
    }
    return("")
  }
  
if(want.trend) {
  # Check for use of internal names in trend
  cc <- check.clashes(internal.names, termsinformula(trend),
                      "the model formula")
  if(cc != "") stop(cc)
  # Default explanatory variables for trend
  glmdata <- data.frame(glmdata, x=P$x, y=P$y)
  if(!is.null(.mpl$MARKS))
    glmdata <- data.frame(glmdata, marks=.mpl$MARKS)
  # 
  if(!is.null(covariates)) {
#   Check for duplication of reserved names
    cc <- check.clashes(reserved.names, names(covariates), "\'covariates\'")
    if(cc != "") stop(cc)
#   Append `covariates.df' to `glmdata'
    glmdata <- data.frame(glmdata,covariates.df)
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

  # List of interaction variable names
  Vnames <- dimnames(V)[[2]]
  
  #   Check for name clashes between the interaction variables
  #   and the formula
  cc <- check.clashes(Vnames, termsinformula(trend), "model formula")
  if(cc != "") stop(cc)
  #   and with the variables in 'covariates'
  if(!is.null(covariates)) {
    cc <- check.clashes(Vnames, names(covariates), "\'covariates\'")
    if(cc != "") stop(cc)
  }

  # OK. append variables.
  glmdata <- data.frame(glmdata, V)   

# Keep only those quadrature points for which the
# conditional intensity is nonzero. 

#KEEP  <- apply(V != -Inf, 1, all)
.mpl$KEEP  <- matrowall(V != -Inf)

.mpl$SUBSET <- .mpl$SUBSET & .mpl$KEEP

if(any(.mpl$Z & !.mpl$KEEP)) {
        howmany <- sum(.mpl$Z & !.mpl$KEEP)
	warning(paste(howmany, "data point(s) are illegal (zero conditional intensity under the model)"))
#	browser()
}

}

##################   D a t a    f r a m e   ###################

# Determine the domain of integration for the pseudolikelihood.

if(correction == "border") {
	bd <- bdist.points(P)
	.mpl$DOMAIN <- (bd >= rbord)
	.mpl$SUBSET <- .mpl$DOMAIN & .mpl$SUBSET
}

glmdata <- data.frame(glmdata, .mpl.SUBSET=.mpl$SUBSET)

#################  F o r m u l a   ##################################

if(!want.trend) trend <- ~1 
trendpart <- paste(as.character(trend), collapse=" ")
rhs <- paste(c(trendpart, Vnames), collapse= "+")
fmla <- paste(".mpl.Y ", rhs)
fmla <- as.formula(fmla)

#### 

return(list(fmla=fmla, glmdata=glmdata, Vnames=Vnames))

}


####################################################################
####################################################################

mpl.get.covariates <- function(covariates, locations, type="") {
  x <- locations$x
  y <- locations$y
  if(is.null(x) || is.null(y)) {
    xy <- xy.coords(locations)
    x <- xy$x
    y <- xy$y
  }
  if(is.null(x) || is.null(y))
    stop("Can't interpret \`locations\' as x,y coordinates")
  n <- length(x)
  if(is.data.frame(covariates)) {
    if(nrow(covariates) != n)
      stop(paste("Number of rows in \`covariates\' != number of", type))
    return(covariates)
  } else if(is.list(covariates)) {
    if(!all(unlist(lapply(covariates, is.im))))
      stop("Some entries in the list \`covariates\' are not images")
    if(any(names(covariates) == ""))
      stop("Some entries in the list \`covariates\' are un-named")
    # look up values of each covariate image at the quadrature points
    values <- lapply(covariates, lookup.im, x=x, y=y, naok=TRUE)
    return(as.data.frame(values))
  } else
    stop("\`covariates\' must be either a data frame or a list of images")
}

