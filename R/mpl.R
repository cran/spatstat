#    mpl.R
#
#	$Revision: 5.32 $	$Date: 2006/06/29 00:31:37 $
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
   ppm(Q=Q, trend=trend, interaction=interaction,
       covariates=data, correction=correction, rbord=rbord,
       use.gam=use.gam, method="mpl")
}

"mpl.engine" <- 
function(Q,
         trend = ~1,
	 interaction = NULL,
         ...,
         covariates = NULL,
	 correction="border",
	 rbord = 0,
         use.gam=FALSE,
         forcefit=FALSE,
         callstring="",
         precomputed=NULL,
         savecomputed=FALSE
) {
  if(!is.null(precomputed$Q)) {
    Q <- precomputed$Q
    X <- precomputed$X
    P <- precomputed$U
  } else {
#
# Determine quadrature scheme from argument Q
#
    if(verifyclass(Q, "ppp", fatal = FALSE)) {
      # point pattern - create default quadrature scheme 
      Q <- quadscheme(Q)   
    } else if(verifyclass(Q, "quad", fatal=FALSE)) {
      # user-supplied quadrature scheme - validate it
      validate.quad(Q, fatal=TRUE, repair=FALSE, announce=TRUE)
    } else 
      stop("First argument Q should be a point pattern or a quadrature scheme")
    
#
# Extract data points
    X <- Q$data
#
# Data and dummy points together
    P <- union.quad(Q)

  }
#
#  
  computed <- if(savecomputed) list(X=X, Q=Q, U=P) else NULL
#
#
# Interpret the call
want.trend <- !is.null(trend) && !identical.formulae(trend, ~1)
want.inter <- !is.null(interaction) && !is.null(interaction$family)

the.version <- list(major=1,
                    minor=9,
                    release=4,
                    date="$Date: 2006/06/29 00:31:37 $")

if(use.gam && exists("is.R") && is.R()) 
  require(mgcv)

  
if(!want.trend && !want.inter && !forcefit) {
  # the model is the uniform Poisson process
  # The MPLE (= MLE) can be evaluated directly
  npts <- X$n
  volume <- area.owin(X$window) * markspace.integral(X)
  lambda <- npts/volume
  theta <- list("log(lambda)"=log(lambda))
  maxlogpl <- if(npts == 0) 0 else npts * (log(lambda) - 1)
  rslt <- list(
               method      = "mpl",
               theta       = theta,
               coef        = theta,
               trend       = NULL,
               interaction = NULL,
               Q           = Q,
               maxlogpl    = maxlogpl,
               internal    = list(computed=computed),
               covariates  = NULL,
	       correction  = correction,
               rbord       = rbord,
               version     = the.version,
               problems    = list())
  class(rslt) <- "ppm"
  return(rslt)
}

        
#################  P r e p a r e    D a t a   ######################
        
prep <- mpl.prepare(Q, X, P, trend, interaction,
                    covariates, want.trend, want.inter, correction, rbord,
                    "quadrature points", callstring,
                    precomputed=precomputed, savecomputed=savecomputed)

fmla <- prep$fmla
glmdata <- prep$glmdata
problems <- prep$problems
computed <- append(computed, prep$computed)
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
             internal     = list(glmfit=FIT, glmdata=glmdata, Vnames=Vnames,
               computed=computed),
             covariates   = covariates,
             correction   = correction,
             rbord        = rbord,
             version      = the.version,
             problems     = problems)
class(rslt) <- "ppm"
return(rslt)
}  


##########################################################################
### /////////////////////////////////////////////////////////////////////
##########################################################################


mpl.prepare <- function(Q, X, P, trend, interaction, covariates, 
                        want.trend, want.inter, correction, rbord,
                        Pname="quadrature points", callstring="",
                        ..., 
                        precomputed=NULL, savecomputed=FALSE) {

  if(missing(want.trend))
    want.trend <- !is.null(trend) && !identical.formulae(trend, ~1)
  if(missing(want.inter))
    want.inter <- !is.null(interaction) && !is.null(interaction$family)
    
  computed <- list()
  problems <- list()
  
  names.precomputed <- names(precomputed)
  

################ C o m p u t e     d a t a  ####################

# Extract covariate values
  if(want.trend && !is.null(covariates)) {
    if("covariates.df" %in% names.precomputed)
      covariates.df <- precomputed$covariates.df
    else 
      covariates.df <- mpl.get.covariates(covariates, P, Pname)
    if(savecomputed)
      computed$covariates.df <- covariates.df
  }
        
### Form the weights and the ``response variable''.

  if("dotmplbase" %in% names.precomputed) 
    .mpl <- precomputed$dotmplbase
  else {
    .mpl <- list()
    .mpl$W <- w.quad(Q)
    .mpl$Z <- is.data(Q)
    .mpl$Y <- .mpl$Z/.mpl$W
    .mpl$MARKS <- marks.quad(Q)  # is NULL for unmarked patterns
    n <- n.quad(Q)
    .mpl$SUBSET <- rep(TRUE, n)
	
    zeroes <- attr(.mpl$W, "zeroes")
    if(!is.null(zeroes))
      .mpl$SUBSET <-  !zeroes
  }

  if(savecomputed)
    computed$dotmplbase <- .mpl
  
  glmdata <- data.frame(.mpl.W = .mpl$W,
                        .mpl.Y = .mpl$Y)
        
    
####################### T r e n d ##############################

  internal.names <- c(".mpl.W", ".mpl.Y", ".mpl.Z", ".mpl.SUBSET",
                        "SUBSET", ".mpl")

  reserved.names <- c("x", "y", "marks", internal.names)

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
  trendvariables <- variablesinformula(trend)
  # Check for use of internal names in trend
  cc <- check.clashes(internal.names, trendvariables, "the model formula")
  if(cc != "") stop(cc)
  # Standard variables
  if("x" %in% trendvariables)
    glmdata <- data.frame(glmdata, x=P$x)
  if("y" %in% trendvariables)
    glmdata <- data.frame(glmdata, y=P$y)
  if(("marks" %in% trendvariables) || !is.null(.mpl$MARKS))
    glmdata <- data.frame(glmdata, marks=.mpl$MARKS)
  #
  # Check covariates
  if(!is.null(covariates)) {
#   Check for duplication of reserved names
    cc <- check.clashes(reserved.names, names(covariates), "\'covariates\'")
    if(cc != "") stop(cc)
#   Take only those covariates that are named in the trend formula
    needed <- names(covariates.df) %in% trendvariables
    covariates.needed <- covariates.df[, needed, drop=FALSE]
#   Append to `glmdata'
    glmdata <- data.frame(glmdata,covariates.needed)
#   Ignore any quadrature points that have NA's in the covariates
    nbg <- is.na(covariates.needed)
    if(any(nbg)) {
      offending <- matcolany(nbg)
      covnames.na <- names(covariates.needed)[offending]
      quadpoints.na <- matrowany(nbg)
      n.na <- sum(quadpoints.na)
      n.tot <- length(quadpoints.na)
      complaint <- paste("Values of the",
                         ngettext(length(covnames.na),
                                  "covariate", "covariates"),
                         paste(sQuote(covnames.na), collapse=", "),
                         "were NA or undefined at",
                         paste(ceiling(100 * n.na/n.tot),
                               "% (", n.na, " out of ", n.tot, ")",
                               sep=""),
                         "of the", Pname)
      warning(paste(complaint,
                    ". Occurred while executing: ",
                    callstring, sep=""),
              call. = FALSE)
      .mpl$SUBSET <-  .mpl$SUBSET & !quadpoints.na
      details <- list(covnames.na   = covnames.na,
                      quadpoints.na = quadpoints.na,
                      print         = complaint)
      problems <- append(problems,
                         list(na.covariates=details))
    }
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

  evaluate <- interaction$family$eval
  if("precomputed" %in% names(formals(evaluate))) {
    # version 1.9-3 onward
    V <- evaluate(X, P, E,
                  interaction$pot,
                  interaction$par,
                  correction, ...,
                  precomputed=precomputed,
                  savecomputed=savecomputed)
    # extract intermediate computation results 
    if(savecomputed)
      computed <- append(computed, attr(V, "computed"))
  } else {
    # Object created by earlier version of ppm.
    # Cannot use precomputed data
    V <- evaluate(X, P, E,
                  interaction$pot,
                  interaction$par,
                  correction)
  }

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
  
  # Check the names are valid as column names in a dataframe
  okVnames <- make.names(Vnames)
  if(any(Vnames != okVnames)) {
    warning("internal error: names of interaction terms contained illegal characters; names have been repaired.")
    Vnames <- okVnames
  }
    
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
        complaint <- paste(howmany,
                           "data point(s) are illegal",
                           "(zero conditional intensity under the model)")
        details <- list(illegal=howmany,
                        print=complaint)
        problems <- append(problems, list(zerolikelihood=details))
        warning(paste(complaint,
                      ". Occurred while executing: ",
                      callstring, sep=""),
                call. = FALSE)
#	browser()
}

}

##################   D a t a    f r a m e   ###################

# Determine the domain of integration for the pseudolikelihood.

if(correction == "border") {

  if("bdP" %in% names.precomputed)
    bd <- precomputed$bdP
  else
    bd <- bdist.points(P)

  if(savecomputed)
    computed$bdP <- bd
  
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

return(list(fmla=fmla, glmdata=glmdata, Vnames=Vnames, problems=problems,
            computed=computed))

}


####################################################################
####################################################################

mpl.get.covariates <- function(covariates, locations, type="locations") {
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

  
  
