#
#
#   rmhcontrol.R
#
#   $Revision: 1.4 $  $Date: 2005/03/10 17:29:10 $
#
#

rmhcontrol <- function(control, ...) {
  UseMethod("rmhcontrol")
}

rmhcontrol.rmhcontrol <- function(control, ...) {
  return(control)
}

rmhcontrol.list <- function(control, ...) {
  argnames <- c("p","q","nrep","expand","periodic",
                "ptypes","fixall", "nverb")
  ok <- argnames %in% names(control) 
  co <- do.call("rmhcontrol.default", control[argnames[ok]])
  return(co)
}

rmhcontrol.default <- function(control=NULL, ..., p=0.9, q=0.5, nrep=5e5,
                        expand=NULL, periodic=FALSE, ptypes=NULL,
                        fixall=FALSE, nverb=0)
{
  if(!is.null(control) || length(list(...)) > 0)
    stop(paste("Syntax should be rmhcontrol(",
               "p, q, nrep, expand, periodic, ptypes, fixall, nverb",
               ")\n with arguments given by name if present"))
  
  # validate arguments
  if(!is.numeric(p) || length(p) != 1
     || p < 0 || p > 1)
    stop("p should be a number in [0,1]")
  if(!is.numeric(q) || length(q) != 1
     || q < 0 || q > 1)
    stop("q should be a number in [0,1]")
  if(!is.numeric(nrep) || length(nrep) != 1
     || nrep <= 1)
    stop("nrep should be an integer > 1")
  if(!is.numeric(nverb) || length(nverb) != 1
     || nverb < 0 || nverb > nrep)
    stop("nverb should be an integer <= nrep")
  if(!is.logical(fixall) || length(fixall) != 1)
    stop("fixall should be a logical value")
  if(!is.logical(periodic) || length(periodic) != 1)
    stop("\`periodic\' should be a logical value")
  if(!is.null(expand) && !(is.numeric(expand) || is.owin(expand)))
    stop("\`expand\' should be either a number, a window, or NULL")


#################################################################
# Conditioning on the number of points?
#  
# cond = 1 <--> no conditioning
# cond = 2 <--> conditioning on n = number of points
# cond = 3 <--> conditioning on the number of points of each type.

  cond    <- 2 - (p<1) + fixall - fixall*(p<1)
  conditioning <- switch(cond, "none", "n.total", "n.each.type")
  
# Warn about silly combination
  if(fixall && p < 1)
	warning("fixall = TRUE conflicts with p < 1. Ignored.\n")

###############################################################  
# `expand' determines expansion of the simulation window

  if(is.null(expand)) 
    force.exp <- force.noexp <- FALSE
  else if(is.owin(expand)) {
    force.exp <- TRUE
    force.noexp <- FALSE
  } else if(is.numeric(expand)) {
    if(expand < 1) {
      warning("parameter \'expand\' < 1 reset to 1")
      expand <- 1
    }
    force.exp <- (expand > 1)
    force.noexp <- !force.exp
  } else stop(paste("Component expand of argument control must",
                    "be either a number, a window, or NULL.\n"))

# No expansion is permitted if we are conditioning on the
# number of points
  
  if(conditioning != "none") {
    if(force.exp)
      stop(paste("When conditioning on the number of points,",
                 "no expansion may be done.\n"))
    else if(is.null(expand)) expand <- 1
  }

# At this stage if there was a reason not to expand the window,
# then expand has been specified (and set to 1).  So if expand has
# not been specified we let it default to the numeric value 2.
# However force.exp = force.noexp = FALSE
# indicating that this was not the user's idea  
  if(is.null(expand)) expand <- 2

###################################################################
# return augmented list  
  out <- list(p=p, q=q, 
              nrep=nrep, nverb=nverb,
              expand=expand, 
              periodic=periodic, 
              ptypes=ptypes,
              fixall=fixall,
              force.exp=force.exp,
              force.noexp=force.noexp,
              cond=cond,
              conditioning=conditioning)
  class(out) <- c("rmhcontrol", class(out))
  return(out)
}

print.rmhcontrol <- function(x, ...) {
  verifyclass(x, "rmhcontrol")

  cat("Metropolis-Hastings algorithm control parameters\n")
  cat(paste("Probability of shift proposal: p =", x$p, "\n"))
  cat(paste("Conditional probability of death proposal: q =", x$q, "\n"))
  switch(x$cond,
         {},
         cat("The total number of points is fixed\n"),
         cat("The number of points of each type is fixed\n"))
  if(!is.null(x$ptypes)) {
    cat("Birth proposal probabilities for each type of point:\n")
    print(x$ptypes)
  }

  cat(paste("Number of M-H iterations: nrep =", x$nrep, "\n"))
  if(x$nverb > 0)
    cat(paste("Progress report every nverb=", x$nverb, "iterations\n"))
  else
    cat("No progress reports (nverb = 0).\n")

  cat("Expand the simulation window? ")
  if(x$force.noexp)
    cat("No.\n")
  else {
    if(x$force.exp)
      cat("Yes:\n")
    else
      cat("Not determined. Default is:\n")
    
    if(is.numeric(x$expand))
      cat(paste("\tarea expansion factor", x$expand, "\n"))
    else 
      print(x$expand)
  }
}
