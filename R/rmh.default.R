#
# $Id: rmh.default.R,v 1.39 2006/05/03 09:14:39 adrian Exp adrian $
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

  if(is.null(start$seed))
    start$seed <- rmhseed()
  
#### Multitype models
  
# Decide whether the model is multitype; if so, find the types.

  types <- rmhResolveTypes(model, start, control)
  ntypes <- length(types)
  mtype <- (ntypes > 1)
  model$types <- types

# If the model is multitype, check that the model parameters agree with types
# and digest them
  
  if(mtype && !is.null(model$check)) 
    model$fortran.par <- model$check(model$par, types)

  
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
      stop("Internal error: unrecognised value of \`expand\'")

  expanded <- (!is.numeric(expand) || (expand > 1))

## Check the fine print   

  if(expanded) {

    if(control$conditioning != "none")
      stop(paste("If we're conditioning on the number of points,",
                 "we cannot clip the result to another window.\n"))

    if(!is.subset.owin(w.clip, w.sim))
      stop("Expanded simulation window does not contain clipping window")
  }


######### Store decisions

  Model <- model
  Start <- start
  Control <- control

  Model$w <- w.clip
  Control$expand <- if(expanded) w.sim else 1
  Control$force.exp <- expanded
  Control$force.noexp <- !expanded

  InfoList <- list(model=Model, start=Start, control=Control)

#########  Make an empty pattern, to be returned in some cases

  empty <- ppp(numeric(0), numeric(0), window=w.clip)
  if(mtype) {
    vide <- factor(types[integer(0)], levels=types)
    empty <- empty %mark% vide
  }
  attr(empty, "info") <- InfoList
  
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
        Control$expand <- w.sim
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

# If we're simulating the empty pattern, return it now
if(npts == 0 && control$conditioning != "none") {
  if(verbose) 
    cat("Initial pattern has 0 points, and simulation is conditional on the number of points - returning an empty pattern\n")
  return(empty)
}

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

# parameter passed to Fortran:  
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
                       x = table(x.start$marks)/x.start$n
                       )
    } else {
      # Validate ptypes
      if(length(ptypes) != ntypes | sum(ptypes) != 1)
        stop("Argument ptypes is mis-specified.\n")
    }
  } 

#==+===+===+===+===+===+===+===+===+===+===+===+===+===+===+===+===+===+===
#
#     S I M U L A T I O N     
#
#==+===+===+===+===+===+===+===+===+===+===+===+===+===+===+===+===+===+===

####### a.s. empty pattern      ###############################
#
#  If beta = 0, the process is almost surely empty
#  
  
  if(all(model$par[["beta"]] < .Machine$double.eps)) {
    if(control$conditioning == "none") {
      # return empty pattern
      if(verbose)
        cat("beta = 0 implies an empty pattern\n")
      return(empty)
    } else 
      stop("beta = 0 implies an empty pattern, but we are simulating conditional on a nonzero number of points")
  }
    
  
#######  Poisson case  ################################
#
#
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
                        x = as.integer(table(x.start$marks)),
  stop("No starting state given; can't condition on fixed number of points"))
               rmpoint(npts.each, intensity, win=w.sim, types=types,
                       verbose=verbose)
             },
             stop("Internal error: control$conditioning unrecognised")
             )
    Xclip <- Xsim[, w.clip]
    attr(Xclip, "info") <- InfoList
    return(Xclip)
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

########################################################################  
#      S t a r t   s i m u l a t i o n
########################################################################

  if(verbose)
    cat("Starting simulation.\nInitial state...")
  
####  Random number initialisation
  
  set.seed(start$seed$build.seed)


#### Build starting state


#### First the marks, if any.
#### The marks must be integers 1 to ntypes, for passing to Fortran
  
  marks <-
    if(!mtype)
      0
    else
      switch(start$given,
             n = {
               # n.start given
               if(control$conditioning=="n.each.type")
                 rep(1:ntypes,n.start)
               else
                 sample(1:ntypes,npts,TRUE,ptypes)
             },
             x = {
               # x.start given
               as.integer(x.start$marks)
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
                      factor(marks,levels=1:ntypes), w.sim, ...)
           x <- xy$x
           y <- xy$y
         })


#######################################################################
#  Start to call Fortran
######################################################################    

# Get the parameters in Fortran-ese
    
  par <- model$fortran.par
  nmbr <- model$fortran.id

# Absorb the constants or vectors `iota' and 'ptypes' into the beta parameters
# This assumes that the beta parameters are the first 'ntypes' entries
# of the parameter vector passed to Fortran.

  par[1:ntypes] <- (iota/ptypes) * par[1:ntypes]
  
# Algorithm control parameters

  nrep <- control$nrep
  cond <- control$cond
  
# If we are simulating a Geyer saturation process we need to set up some
# ``auxiliary information''.

  need.aux <- model$need.aux
  aux <-
    if(!need.aux)
      0
    else
      .Fortran(
               "initaux",
               nmbr=as.integer(nmbr),
               par=as.double(par),
               period=as.double(period),
               x=as.double(x),
               y=as.double(y),
               npts=as.integer(npts),
               aux=integer(npts),
               PACKAGE="spatstat"
               )$aux

# The vectors x and y (and perhaps marks) which hold the generated
# process may grow.  We need to allow storage space for them to grow
# in.  Unless we are conditioning on the number of points, we have no
# real idea how big they will grow.  Hence we start off with storage
# space which has at least twice the length of the ``initial state'', and
# structure things so that the storage space may be incremented
# without losing the ``state'' which has already been generated.

  if(cond == 1) {
    nincr <- npad <- max(npts, 50)
    padding <- numeric(npad)
    x <- c(x, padding)
    y <- c(y, padding)
    if(ntypes>1) marks <- c(marks, padding)
    if(need.aux) aux   <- c(aux,   padding)
  } else {
    nincr <- npad <- 0
  }
  npmax <- npts + npad
  mrep  <- 1

  if(verbose)
    cat("Proposal points...")
           
# If the pattern is multitype, generate the mark proposals.
  mprop <- if(ntypes>1)
    sample(1:ntypes,nrep,TRUE,prob=ptypes) else 0
		
# Generate the ``proposal points'' in the expanded window.
  xy <-
    if(trendy)
      rpoint.multi(nrep,trend,tmax,factor(mprop, levels=1:ntypes),w.sim,...)
    else
      runifpoint(nrep,w.sim)
  xprop <- xy$x
  yprop <- xy$y

  if(verbose)
    cat("Start simulation.\n")

# Determine Fortran subroutine name (this is not very safe ...)
  mhname <- paste("mh", nmbr, sep="")
  
# The repetition is to allow the storage space to be incremented if
# necessary.
  repeat {
# Call the Metropolis-Hastings simulator:
    rslt <- .Fortran(
                     mhname,
                     par=as.double(par),
                     period=as.double(period),
                     xprop=as.double(xprop),
                     yprop=as.double(yprop),
                     mprop=as.integer(mprop),
                     ntypes=as.integer(ntypes),
                     iseed=as.integer(start$seed$iseed),
                     nrep=as.integer(nrep),
                     mrep=as.integer(mrep),
                     p=as.double(control$p),
                     q=as.double(control$q),
                     npmax=as.integer(npmax),
                     nverb=as.integer(control$nverb),
                     x=as.double(x),
                     y=as.double(y),
                     marks=as.integer(marks),
                     aux=as.integer(aux),
                     npts=as.integer(npts),
                     fixall=as.logical(control$fixall),
                     PACKAGE="spatstat"
                     )

    npts <- rslt$npts
    mrep <- rslt$mrep
    
# If mrep > nrep we have completed the nrep Metropolis Hastings steps.
# Exit the loop.
    
    if(mrep > nrep) break

# If not, then we came back early because we were about to run out
# of storage space.  Increase the storage space and re-call the
# methas subroutine.

    # Internal consistency check
    if(npts < npmax)
	stop(paste("Internal error; Fortran code exited unexpectedly\n",
                   "(without completing nrep steps and without reaching\n",
                   "the limit of storage capacity)\n"))
    
    if(verbose) {
      cat('Number of points equal to ',npmax,';\n',sep='')
      cat('increasing storage space and continuing.\n')
    }
    npmax <- npmax + nincr
    x     <- c(rslt$x,numeric(nincr))
    y     <- c(rslt$y,numeric(nincr))
    marks <- if(ntypes>1) c(rslt$marks,numeric(nincr)) else 0
    aux   <- if(need.aux) c(rslt$aux,numeric(nincr)) else 0
    iseed <- rslt$iseed
  }

  ###### END OF LOOP ###################
  
  # Extract the point pattern returned from Fortran
  
  npts <- rslt$npts
  indices <- if(npts == 0) numeric(0) else (1:npts)
  x <- rslt$x[indices]
  y <- rslt$y[indices]
  xxx <- ppp(x=x, y=y, window=w.sim)
  if(mtype) {
    marx <- factor(rslt$marks[indices],levels=types)
    xxx <- xxx %mark% marx
  } 

# Now clip the pattern to the ``clipping'' window:
  xxx <- xxx[,w.clip]

# Append to the result information about how it was generated.
  attr(xxx, "info") <- InfoList
  return(xxx)
}


