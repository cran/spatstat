#
#
#   rmhResolveTypes.R
#
#   $Revision: 1.1 $   $Date: 2005/03/05 01:11:00 $
#
#
rmhResolveTypes <- function(model, start, control) {

# Decide whether a multitype point process is to be simulated.
# If so, determine the vector of types.

  verifyclass(model, "rmhmodel")
  verifyclass(start, "rmhstart")
  verifyclass(control, "rmhcontrol")

# Different ways of specifying types

  types.model <- model$types
  types.start <- if(start$given=="x" && is.marked(x.start <- start$x.start))
                     levels(x.start$marks) else NULL
  
# Check for inconsistencies  
  if(!is.null(types.model) && !is.null(types.start))
    if(!identical(all.equal(types.model, types.start), TRUE))
      stop("marks in control$x.start do not match model$types")
  
  types.given <- if(!is.null(types.model)) types.model else types.start
  types.given.source <-
    if(!is.null(types.model)) "model$types" else "marks of x.start"
  
# Different ways of specifying/implying the number of types
  
  ntypes.beta <- length(model$par$beta)
  ntypes.ptypes <- length(control$ptypes)
  ntypes.nstart <- if(start$given == "n") length(start$n.start) else 0
  ntypes.trend <- if(is.list(model$trend)) length(model$trend) else 1
  
# Check for inconsistencies (only for numbers > 1)

  nty <- c(ntypes.beta, ntypes.ptypes, ntypes.nstart)
  nam <- c("model$par$beta", "control$ptypes", "start$n.start", "model$trend")
  give <- (nty > 1)
  if(!any(give))
    ntypes.given <- 1
  else {
    if(length(unique(nty[give])) > 1)
      stop(paste("Mismatch in lengths of",
               paste(nam[give], collapse=", ")))
    ntypes.given <- unique(nty[give])
    ntypes.given.source <- (nam[give])[1]
  } 

# Check types.given and ntypes.given 

  if(!is.null(types.given) && ntypes.given > 1)
    if(length(types.given) != ntypes.given)
      stop(paste("Mismatch between number of types in",
                 types.given.source,
                 "and length of",
                 ntypes.given.source))

# Finally determine the types
  
  if(model$multitype.interact) {
    # There MUST be a types vector
    types <- if(!is.null(types.given)) types.given
             else if(ntypes.given > 1) 1:ntypes.given
             else stop("Cannot determine types for multitype process")
  } else {
    types <- if(!is.null(types.given)) types.given
             else if(ntypes.given > 1) 1:ntypes.given
             else 1
  }

  ntypes <- length(types)
  
# If we are conditioning on the number of points of each type,
# make sure the starting state is appropriate

  if(control$conditioning == "n.each.type") {
    if(start$given == "n" && ntypes.nstart != ntypes)
      stop("Length of control$n.start not equal to number of types.\n")
    else if(start$given == "x" && ntypes.start != ntypes) 
      stop("Marks of control$x.start do not match number of types.\n")
  }
  
# Warn about a silly value of fixall:
  if(control$fixall & ntypes==1)
	warning("fixall = TRUE conflicts with ntypes = 1. Ignored. \n")

  return(types)
}

  
