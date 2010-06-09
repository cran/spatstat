#
#
#   rmhcontrol.R
#
#   $Revision: 1.10 $  $Date: 2010/06/05 09:49:46 $
#
#

rmhcontrol <- function(...) {
  UseMethod("rmhcontrol")
}

rmhcontrol.rmhcontrol <- function(...) {
  argz <- list(...)
  if(length(argz) == 1)
    return(argz[[1]])
  stop("Arguments not understood")
}

rmhcontrol.list <- function(...) {
  argz <- list(...)
  nama <- names(argz)
  if(length(argz) == 1 && !any(nzchar(nama)))
    do.call("rmhcontrol.default", argz[[1]])
  else
    do.call.matched("rmhcontrol.default", argz)
}

rmhcontrol.default <- function(..., p=0.9, q=0.5, nrep=5e5,
                        expand=NULL, periodic=FALSE, ptypes=NULL,
                        x.cond=NULL, fixall=FALSE, nverb=0)
{
  argh <- list(...)
  nargh <- length(argh)
  if(nargh > 0) {
    # allow rmhcontrol(NULL), otherwise flag an error
    if(!(nargh == 1 && is.null(argh[[1]])))
      stop(paste("Unrecognised arguments; syntax should be rmhcontrol(",
                 "p, q, nrep, expand, periodic, ptypes, x.cond, fixall, nverb",
                 ") with arguments given explicitly by name"))
  }
  # validate arguments
  if(!is.numeric(p) || length(p) != 1
     || p < 0 || p > 1)
    stop("p should be a number in [0,1]")
  if(!is.numeric(q) || length(q) != 1
     || q < 0 || q > 1)
    stop("q should be a number in [0,1]")
  if(!is.numeric(nrep) || length(nrep) != 1
     || nrep < 1)
    stop("nrep should be an integer >= 1")
  if(!is.numeric(nverb) || length(nverb) != 1
     || nverb < 0 || nverb > nrep)
    stop("nverb should be an integer <= nrep")
  if(!is.logical(fixall) || length(fixall) != 1)
    stop("fixall should be a logical value")
  if(!is.logical(periodic) || length(periodic) != 1)
    stop(paste(sQuote("periodic"), "should be a logical value"))
  if(!is.null(expand) && !(is.numeric(expand) || is.owin(expand)))
    stop(paste(sQuote("expand"),
               "should be either a number, a window, or NULL"))

#################################################################
# Conditioning on point configuration
#
# condtype = "none": no conditioning
# condtype = "Palm": conditioning on the presence of specified points
# condtype = "window": conditioning on the configuration in a subwindow
#
  if(is.null(x.cond)) {
    condtype <- "none"
    n.cond <- NULL
  } else if(is.ppp(x.cond)) {
    condtype <- "window"
    n.cond <- x.cond$n
  } else if(is.data.frame(x.cond)) {
    if(ncol(x.cond) %in% c(2,3)) {
      condtype <- "Palm"
      n.cond <- nrow(x.cond)
    } else stop("Wrong number of columns in data frame x.cond")
  } else if(is.list(x.cond)) {
    if(length(x.cond) %in% c(2,3)) {
      x.cond <- as.data.frame(x.cond)
      condtype <- "Palm"
      n.cond <- nrow(x.cond)
    } else stop("Wrong number of components in list x.cond")
  } else stop("Unrecognised format for x.cond")

  if(condtype == "Palm" && n.cond == 0) {
    warning(paste("Ignored empty configuration x.cond;",
                  "conditional (Palm) simulation given an empty point pattern",
                  "is equivalent to unconditional simulation"))
    condtype <- "none"
    x.cond <- NULL
    n.cond <- NULL
  }
    
#################################################################
# Fixing the number of points?
#  
# fixcode = 1 <--> no conditioning
# fixcode = 2 <--> conditioning on n = number of points
# fixcode = 3 <--> conditioning on the number of points of each type.

  fixcode    <- 2 - (p<1) + fixall - fixall*(p<1)
  fixing <- switch(fixcode, "none", "n.total", "n.each.type")
  
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
      warning(paste("parameter", sQuote("expand"), " < 1; reset to 1"))
      expand <- 1
    }
    force.exp <- (expand > 1)
    force.noexp <- !force.exp
  } else stop(paste("Component expand of argument control must",
                    "be either a number, a window, or NULL.\n"))

# No expansion is permitted if we are conditioning on the
# number of points
  
  if(fixing != "none") {
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
  if(is.null(expand))
    expand <- spatstat.options("expand")

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
              fixcode=fixcode,
              fixing=fixing,
              condtype=condtype,
              x.cond=x.cond)
  class(out) <- c("rmhcontrol", class(out))
  return(out)
}

print.rmhcontrol <- function(x, ...) {
  verifyclass(x, "rmhcontrol")

  cat("Metropolis-Hastings algorithm control parameters\n")
  cat(paste("Probability of shift proposal: p =", x$p, "\n"))
  if(x$fixing == "none") {
    cat(paste("Conditional probability of death proposal: q =", x$q, "\n"))
    if(!is.null(x$ptypes)) {
      cat("Birth proposal probabilities for each type of point:\n")
      print(x$ptypes)
    }
  }
  switch(x$fixing,
         none={},
         n.total=cat("The total number of points is fixed\n"),
         n.each.type=cat("The number of points of each type is fixed\n"))
  switch(x$condtype,
         none={},
         window={
           cat(paste("Conditional simulation given the",
                     "configuration in a subwindow\n"))
           print(x$x.cond$window)
         },
         Palm={
           cat("Conditional simulation of Palm type\n")
         })
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
