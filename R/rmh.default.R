#
# $Id: rmh.default.R,v 1.58 2009/10/22 18:42:20 adrian Exp adrian $
#
rmh.default <- function(model,start=NULL,control=NULL, verbose=TRUE, ...) {
#
# Function rmh.  To simulate realizations of 2-dimensional point
# patterns, given the conditional intensity function of the 
# underlying process, via the Metropolis-Hastings algorithm.
#
#==+===+===+===+===+===+===+===+===+===+===+===+===+===+===+===+===+===+===
#
#     V A L I D A T E
#
#==+===+===+===+===+===+===+===+===+===+===+===+===+===+===+===+===+===+===
  
  if(verbose)
    cat("Checking arguments..")
  
# validate arguments and fill in the defaults
  
  model <- rmhmodel(model)
  start <- rmhstart(start)
  control <- rmhcontrol(control)

#### Multitype models
  
# Decide whether the model is multitype; if so, find the types.

  types <- rmhResolveTypes(model, start, control)
  ntypes <- length(types)
  mtype <- (ntypes > 1)
  model$types <- types

# If the model is multitype, check that the model parameters agree with types
# and digest them
  
  if(mtype && !is.null(model$check)) 
    model$C.par <- model$check(model$par, types)

  
######## Check for illegal combinations of model, start and control  ########

  # No expansion can be done if we are using x.start

  if(start$given == "x") {
    if(control$force.exp)
      stop("Cannot expand window when using x.start.\n")
    else 
      control$expand <- 1
  }

# Warn about a silly value of fixall:
  if(control$fixall & ntypes==1) {
    warning("control$fixall applies only to multitype processes. Ignored. \n")
    control$fixall <- FALSE
    if(control$conditioning == "n.each.type")
      control$conditioning <- "n.total"
  }

  
#==+===+===+===+===+===+===+===+===+===+===+===+===+===+===+===+===+===+===
#
#     M O D E L   P A R A M E T E R S
#
#==+===+===+===+===+===+===+===+===+===+===+===+===+===+===+===+===+===+===

#######  Determine windows  ################################

  if(verbose)
    cat("determining simulation windows...")
  
# these may be NULL  
  w.model <- model$w
  x.start <- start$x.start
  trend <- model$trend

  trendy <- !is.null(trend)
  
# window implied by trend image, if any
  
  w.trend <- 
    if(is.im(trend))
      as.owin(trend)
    else if(is.list(trend) && any(ok <- unlist(lapply(trend, is.im))))
      as.owin((trend[ok])[[1]])
    else NULL
 
##  Clipping window (for final result)
  
  w.clip <-
    if(!is.null(model$w))
      model$w
    else if(control$expand == 1) {
      if(start$given == "x" && is.ppp(x.start))
        x.start$window
      else if(is.owin(w.trend))
        w.trend
    } else NULL

  if(!is.owin(w.clip))
    stop("Unable to determine window for pattern")

  
##  Simulation window

  expand <- control$expand

  w.sim <-
    if(is.owin(expand))
      expand
    else if(is.numeric(expand))
      expand.owin(w.clip, expand)
    else
      stop(paste("Internal error: unrecognised value of",
                 sQuote("expand")))

  expanded <- (!is.numeric(expand) || (expand > 1))

## Check the fine print   

  if(expanded) {

    if(control$conditioning != "none")
      stop(paste("If we're conditioning on the number of points,",
                 "we cannot clip the result to another window.\n"))

    if(!is.subset.owin(w.clip, w.sim))
      stop("Expanded simulation window does not contain clipping window")
  }


#######  Trend  ################################
  
# Check that the expanded window fits inside the window
# upon which the trend(s) live if there are trends and
# if any trend is given by an image.

  if(expanded && !is.null(trend)) {
    trends <- if(!is.list(trend)) list(trend) else trend
    images <- unlist(lapply(trends, is.im))
    if(any(images)) {
      iwindows <- lapply(trends[images], as.owin)
      nimages <- length(iwindows)
      misfit <- !unlist(lapply(iwindows,
                               function(x,w) { is.subset.owin(w,x) },
                               w = w.sim))
      nmisfit <- sum(misfit)
      if(nmisfit > 1) 
        stop(paste("Expanded simulation window is not contained in",
                   "several of the trend windows.\n",
                   "Bailing out.\n"))
      else if(nmisfit == 1) {
        warning(paste("Expanded simulation window is not contained in",
                      if(nimages == 1) "the trend window.\n"
                      else "one of the trend windows.\n",
                      "Expanding to this trend window (only).\n"))
        w.sim <- iwindows[misfit]
      }
    }
  }


  
#==+===+===+===+===+===+===+===+===+===+===+===+===+===+===+===+===+===+===
#
#     S T A R T I N G      S T A T E
#
#==+===+===+===+===+===+===+===+===+===+===+===+===+===+===+===+===+===+===

###################### Starting state data ############################

  switch(start$given,
         none={
           stop("No starting state given")
         },
         x = {
           # x.start was given
           # coerce it to a ppp object
           if(!is.ppp(x.start))
             x.start <- as.ppp(x.start, w.sim)
           npts <- x.start$n
         },
         n = {
           # n.start was given
           n.start <- start$n.start
           # Adjust the number of points in the starting state in accordance
           # with the expansion that has occurred.  
           if(expanded)
             n.start <- ceiling(n.start * area.owin(w.sim)/area.owin(w.clip))
           #
           npts <- sum(n.start) # The ``sum()'' is redundant if n.start
                                # is scalar; no harm, but.
         },
         stop("Internal error: start$given unrecognized"))


#==+===+===+===+===+===+===+===+===+===+===+===+===+===+===+===+===+===+===
#
#     C O N T R O L    P A R A M E T E R S
#
#==+===+===+===+===+===+===+===+===+===+===+===+===+===+===+===+===+===+===

###################  Periodic boundary conditions #########################
  
# If periodic is TRUE we have to be simulating in a rectangular window.

  periodic <- control$periodic
  
  if(periodic && w.sim$type != "rectangle")
      stop("Need rectangular window for periodic simulation.\n")

# parameter passed to C:  
  period <-
    if(periodic)
      c(diff(w.sim$xrange), diff(w.sim$yrange))
    else
      c(-1,-1)



#### vector of proposal probabilities 

  if(!mtype) 
    ptypes <- 1
  else {
    ptypes <- control$ptypes
    if(is.null(ptypes)) {
      # default values
      ptypes <- switch(start$given,
                       none = ,
                       n = rep(1/ntypes,ntypes),
                       x = table(marks(x.start, dfok=FALSE))/x.start$n
                       )
    } else {
      # Validate ptypes
      if(length(ptypes) != ntypes | sum(ptypes) != 1)
        stop("Argument ptypes is mis-specified.\n")
    }
  } 


  
########################################################################
#  Normalising constant for proposal density
# 
# Integral of trend over the expanded window (or area of window):
# Iota == Integral Of Trend (or) Area.
  
  if(trendy) {
    if(verbose)
      cat("Evaluating trend integral...")
    tlist <- if(is.function(trend) || is.im(trend)) list(trend) else trend
    tsummaries <- lapply(tlist,
                         function(x, w) {
                           tmp  <- as.im(x, w)[w, drop=FALSE]
                           return(summary(tmp))
                         },
                         w=w.sim)
    nbg  <- unlist(lapply(tsummaries, function(x) { x$min < 0 }))
    if(any(nbg))
      stop("Trend has negative values")
    iota <- unlist(lapply(tsummaries, function(x) { x$integral }))
    tmax <- unlist(lapply(tsummaries, function(x) { x$max }))
  } else {
    iota <- area.owin(w.sim)
    tmax <- NULL
  }

  
#==+===+===+===+===+===+===+===+===+===+===+===+===+===+===+===+===+===+===
#
#     A.S. EMPTY PROCESS
#
#==+===+===+===+===+===+===+===+===+===+===+===+===+===+===+===+===+===+===

  a.s.empty <- FALSE
  
#
#  Empty pattern, simulated conditional on n
#  
  if(npts == 0 && control$conditioning != "none") {
    if(verbose) 
      cat("Initial pattern has 0 points, and simulation is conditional on the number of points - returning an empty pattern\n")
    a.s.empty <- TRUE
  } 

#
#  If beta = 0, the process is almost surely empty
#  
  
  if(all(model$par[["beta"]] < .Machine$double.eps)) {
    if(control$conditioning == "none") {
      # return empty pattern
      if(verbose)
        cat("beta = 0 implies an empty pattern\n")
      a.s.empty <- TRUE
    } else 
      stop("beta = 0 implies an empty pattern, but we are simulating conditional on a nonzero number of points")
  }

  if(a.s.empty) {
    # create empty pattern, to be returned
    empty <- ppp(numeric(0), numeric(0), window=w.clip)
    if(mtype) {
      vide <- factor(types[integer(0)], levels=types)
      empty <- empty %mark% vide
    }
  }

#==+===+===+===+===+===+===+===+===+===+===+===+===+===+===+===+===+===+===
#
#     PACKAGE UP
#
#==+===+===+===+===+===+===+===+===+===+===+===+===+===+===+===+===+===+===

######### Store decisions

  Model <- model
  Start <- start
  Control <- control

  Model$w <- w.clip
  Model$types <- types
  
  Control$expand <- if(expanded) w.sim else 1
  Control$force.exp <- expanded
  Control$force.noexp <- !expanded

  Control$internal <- list(w.sim=w.sim,
                           ptypes=ptypes,
                           period=period)

  Model$internal <- list(a.s.empty=a.s.empty,
                         empty=if(a.s.empty) empty else NULL,
                         mtype=mtype,
                         trendy=trendy,
                         iota=iota,
                         tmax=tmax)

  Start$internal <- list(npts=npts)

  InfoList <- list(model=Model, start=Start, control=Control)
  class(InfoList) <- c("rmhInfoList", class(InfoList))

  # go
  rmhEngine(InfoList, verbose=verbose, kitchensink=TRUE, ...)
}


#---------------  rmhEngine -------------------------------------------
#
# This is the interface to the C code.
#
# InfoList is a list of pre-digested, validated arguments
# obtained from rmh.default.
#
# This function is called by rmh.default to generate one simulated
# realisation of the model.
# It's called repeatedly by ho.engine and qqplot.ppm to generate multiple
# realisations (saving time by not repeating the argument checking
# in rmh.default).

# arguments:  
# kitchensink: whether to tack InfoList on to the return value as an attribute
# preponly: whether to just return InfoList without simulating
#
#   rmh.default digests arguments and calls rmhEngine with kitchensink=T
#
#   qqplot.ppm first gets InfoList by calling rmh.default with preponly=T
#              (which digests the model arguments and calls rmhEngine
#               with preponly=T, returning InfoList),
#              then repeatedly calls rmhEngine(InfoList) to simulate.
#
# -------------------------------------------------------

rmhEngine <- function(InfoList, ...,
                       verbose=FALSE, kitchensink=FALSE,
                       preponly=FALSE) {
# Internal Use Only!
# This is the interface to the C code.

  if(!inherits(InfoList, "rmhInfoList"))
    stop("data not in correct format for internal function rmhEngine")

  
  if(preponly)
    return(InfoList)

  model <- InfoList$model
  start <- InfoList$start
  control <- InfoList$control

  w.sim <- control$internal$w.sim
  w.clip <- model$w

  types <- model$types
  ntypes <- length(types)
  
  ptypes <- control$internal$ptypes
  period <- control$internal$period

  mtype <- model$internal$mtype

  trend <- model$trend
  trendy <- model$internal$trendy
  iota <- model$internal$iota
  tmax <- model$internal$tmax

  npts <- start$internal$npts

  n.start <- start$n.start
  x.start <- start$x.start
  
#==+===+===+===+===+===+===+===+===+===+===+===+===+===+===+===+===+===+===
#
#     E M P T Y   P A T T E R N
#
#==+===+===+===+===+===+===+===+===+===+===+===+===+===+===+===+===+===+===

  if(model$internal$a.s.empty) {
    empty <- model$internal$empty
    attr(empty, "info") <- InfoList
    return(empty)
  }
  
#==+===+===+===+===+===+===+===+===+===+===+===+===+===+===+===+===+===+===
#
#     S I M U L A T I O N     
#
#==+===+===+===+===+===+===+===+===+===+===+===+===+===+===+===+===+===+===

#############################################
####  
####  Random number seed: initialisation & capture
####  
#############################################  
  
  if(!exists(".Random.seed"))
    runif(1)

  saved.seed <- .Random.seed
  

#############################################
####  
####  Poisson case
####  
#############################################  
  
  if(model$cif == 'poisson') {
    intensity <- if(!trendy) model$par[["beta"]] else model$trend
    Xsim <-
      switch(control$conditioning,
             none= {
               # Poisson process 
               if(!mtype)
                 rpoispp(intensity, win=w.sim, ...)
               else
                 rmpoispp(intensity, win=w.sim, types=types)
             },
             n.total = {
               # Binomial/multinomial process with fixed total number of points
               if(!mtype) 
                 rpoint(npts, intensity, win=w.sim, verbose=verbose)
               else
                 rmpoint(npts, intensity, win=w.sim, types=types,
                         verbose=verbose)
             },
             n.each.type = {
               # Multinomial process with fixed number of points of each type
               npts.each <-
                 switch(start$given,
                        n = n.start,
                        x = as.integer(table(marks(x.start, dfok=FALSE))),
  stop("No starting state given; can't condition on fixed number of points"))
               rmpoint(npts.each, intensity, win=w.sim, types=types,
                       verbose=verbose)
             },
             stop("Internal error: control$conditioning unrecognised")
             )
    Xclip <- Xsim[w.clip]
    attr(Xclip, "info") <- InfoList
    return(Xclip)
  }

  
########################################################################  
#      M e t r o p o l i s  H a s t i n g s    s i m u l a t i o n
########################################################################

  if(verbose)
    cat("Starting simulation.\nInitial state...")
  

#### Build starting state


#### First the marks, if any.
#### The marks must be integers 0 to (ntypes-1) for passing to C

  Ctypes <- if(mtype) 0:(ntypes-1) else 0
  
  Cmarks <-
    if(!mtype)
      0
    else
      switch(start$given,
             n = {
               # n.start given
               if(control$conditioning=="n.each.type")
                 rep(Ctypes,n.start)
               else
                 sample(Ctypes,npts,TRUE,ptypes)
             },
             x = {
               # x.start given
               as.integer(marks(x.start, dfok=FALSE))-1
             },
             stop("internal error: start$given unrecognised")
             )
#
# First the x, y coordinates
#
  switch(start$given,
         x = {
           x <- x.start$x
           y <- x.start$y
         },
         n = {
           xy <-
             if(!trendy)
               runifpoint(npts, w.sim, ...)
             else
               rpoint.multi(npts, trend, tmax,
                      factor(Cmarks,levels=Ctypes), w.sim, ...)
           x <- xy$x
           y <- xy$y
         })


#######################################################################
#  Set up C call
######################################################################    

# Determine the name of the cif in the C code

  C.id <- model$C.id
  
# Get the parameters in C-ese
    
  par <- model$C.par

# Absorb the constants or vectors `iota' and 'ptypes' into the beta parameters
# This assumes that the beta parameters are the first 'ntypes' entries
# of the parameter vector passed to Fortran.

  par[1:ntypes] <- (iota/ptypes) * par[1:ntypes]
  
# Algorithm control parameters

  p      <- control$p
  q      <- control$q
  nrep   <- control$nrep
  cond   <- control$cond
  fixall <- control$fixall
  nverb  <- control$nverb
  
# The vectors x and y (and perhaps marks) which hold the generated
# process may grow.  We need to allow storage space for them to grow
# in.  Unless we are conditioning on the number of points, we have no
# real idea how big they will grow.  Hence we start off with storage
# space which has at least twice the length of the ``initial state''.
# If more space is needed, it is re-allocated in the C code.

  if(cond == 1) {
    npad <- max(npts, 50)
    padding <- numeric(npad)
    x <- c(x, padding)
    y <- c(y, padding)
    if(mtype) Cmarks <- c(Cmarks, padding)
  } else {
    npad <- 0
  }
  npmax <- npts + npad
  mrep  <- 1

  if(verbose)
    cat("Proposal points...")
           
# If the pattern is multitype, generate the mark proposals (0 to ntypes-1)
  Cmprop <- if(mtype) sample(Ctypes,nrep,TRUE,prob=ptypes) else 0
		
# Generate the ``proposal points'' in the expanded window.
  xy <-
    if(trendy)
      rpoint.multi(nrep,trend,tmax,factor(Cmprop, levels=Ctypes),w.sim,...)
    else
      runifpoint(nrep,w.sim)
  xprop <- xy$x
  yprop <- xy$y

  if(verbose)
    cat("Start simulation.\n")

# Call the Metropolis-Hastings simulator:
  storage.mode(C.id) <- "character"
  storage.mode(par) <- storage.mode(period) <- "double"
  storage.mode(xprop) <- storage.mode(yprop) <- "double"
  storage.mode(Cmprop) <- "integer"
  storage.mode(ntypes) <- "integer"
  storage.mode(nrep) <- "integer"
  storage.mode(p) <- storage.mode(q) <- "double"
  storage.mode(nverb) <- "integer"
  storage.mode(x) <- storage.mode(y) <- "double"
  storage.mode(Cmarks) <- "integer"
  storage.mode(fixall) <- "integer"

  out <- .Call("methas",
               C.id,
               par,
               period,
               xprop, yprop, Cmprop,
               ntypes, nrep,
               p, q,
               nverb,
               x, y,
               Cmarks, fixall,
               PACKAGE="spatstat")
  
  # Extract the point pattern returned from C

  X <- ppp(x=out[[1]], y=out[[2]], window=w.sim, check=FALSE)
  if(mtype) {
    # convert integer marks from C to R
    marx <- factor(out[[3]]+1, levels=1:ntypes)
    # then restore original type levels
    levels(marx) <- types
    # glue to points
    marks(X) <- marx
  }

# Now clip the pattern to the ``clipping'' window:
  if(!control$force.noexp)
    X <- X[w.clip]

# Append to the result information about how it was generated.
  if(kitchensink) {
    attr(X, "info") <- InfoList
    attr(X, "seed") <- saved.seed
  }
  
  return(X)
}


