#
#	breakpts.S
#
#	A simple class definition for the specification
#       of histogram breakpoints in the special form we need them.
#
#	even.breaks()
#
#	$Revision: 1.5 $	$Date: 2001/11/16 09:03:00 $
#
#
#       Other functions in this directory use the standard Splus function
#	hist() to compute histograms of distance values.
#       One argument of hist() is the vector 'breaks'
#	of breakpoints for the histogram cells. 
#
#       The breakpoints must
#            (a) span the range of the data
#            (b) be given in increasing order
#            (c) satisfy breaks[2] = 0,
#
#	The function make.even.breaks() will create suitable breakpoints.
#
#       Condition (c) means that the first histogram cell has
#       *right* endpoint equal to 0.
#
#       Since all our distance values are nonnegative, the effect of (c) is
#       that the first histogram cell counts the distance values which are
#       exactly equal to 0. Hence F(0), the probability P{X = 0},
#       is estimated without a discretisation bias.
#
#	We assume the histograms have followed the default counting rule
#	in hist(), which is such that the k-th entry of the histogram
#	counts the number of data values in 
#		I_k = ( breaks[k],breaks[k+1] ]	for k > 1
#		I_1 = [ breaks[1],breaks[2]   ]
#
#	The implementations of estimators of c.d.f's in this directory
#       produce vectors of length = length(breaks)-1
#       with value[k] = estimate of F(breaks[k+1]),
#       i.e. value[k] is an estimate of the c.d.f. at the RIGHT endpoint
#       of the kth histogram cell.
#
#       An object of class 'breakpts' contains:
#
#              $val     the actual breakpoints
#              $max     the maximum value (= last breakpoint)
#              $ncells  total number of histogram cells
#              $r       right endpoints, r = val[-1]
#              $even    logical = TRUE if cells known to be evenly spaced
#              $npos    number of histogram cells on the positive halfline
#                        = length(val) - 2,
#                       or NULL if cells not evenly spaced
#              $step    histogram cell width
#                       or NULL if cells not evenly spaced
#       
# --------------------------------------------------------------------
breakpts <- function(val, maxi, even=F, npos=NULL, step=NULL) {
  out <- list(val=val, max=maxi, ncells=length(val)-1, r = val[-1],
              even=even, npos=npos, step=step)
  class(out) <- "breakpts"
  out
}

"make.even.breaks" <- 
function(bmax, npos, bstep) {
        if(missing(bstep) && missing(npos))
          stop("Must specify either \'bstep\' or \'npos\'")
        if(!missing(npos)) {
          bstep <- bmax/npos
          val <- seq(0, bmax, length=npos+1)
          val <- c(-bstep,val)
        } else {
          val <- seq(0, bmax, by=bstep)
          val <- c(-bstep,val)
          npos <- length(val) - 2
        }
        breakpts(val, bmax, T, npos, bstep)
}

"as.breakpts" <- function(...) {

  XL <- list(...)

  if(length(XL) == 1) {
    # single argument
    X <- XL[[1]]

    if(!is.null(class(X)) && class(X) == "breakpts")
    # X already in correct form
      return(X)
  
    if(is.vector(X) && length(X) > 2) {
    # it's a vector
      if(X[2] != 0)
        stop("breakpoints do not satisfy breaks[2] = 0")
      steps <- diff(X)
      if(all(steps == mean(steps)))
        # equally spaced
        return(breakpts(X, max(X), T, length(X)-2, steps[1]))
      else
        # unknown spacing
        return(breakpts(X, max(X), F))
    }
  } else {

    # There are multiple arguments.
  
    # exactly two arguments - interpret as even.breaks()
    if(length(XL) == 2)
      return(make.even.breaks(XL[[1]], XL[[2]]))

    # two arguments 'max' and 'npos'
  
    if(!is.null(XL$max) && !is.null(XL$npos))
      return(make.even.breaks(XL$max, XL$npos))

    # otherwise
    stop("Don't know how to convert these data to breakpoints")
  }
  # never reached
}


check.hist.lengths <- function(hist, breaks) {
  verifyclass(breaks, "breakpts")
  nh <- length(hist)
  nb <- breaks$ncells
  if(nh != nb)
    stop(paste("Length of histogram =", nh,
               "not equal to number of histogram cells =", nb))
}

breakpts.from.r <- function(r) {
        if(r[1] != 0)
          stop("First r value must be 0")
        dr <- r[2] - r[1]
        b <- c(-dr, r)
        return(as.breakpts(b))
}

handle.r.b.args <- function(r=NULL, breaks=NULL, window, eps=NULL) {

        if(!is.null(r) && !is.null(breaks))
          stop("Do not specify both \'r\' and \'breaks\'")
  
        if(!is.null(breaks)) {
          breaks <- as.breakpts(breaks)
        } else if(!is.null(r)) {
          breaks <- breakpts.from.r(r)
	} else {
          # both 'r' and 'breaks' are missing
          if(is.null(eps)) {
            if(!is.null(window$xstep))
              eps <- window$xstep
            else 
              eps <- diff(window$xrange)/100
          }
          # warning(paste("step size for argument \'r\' defaults to", eps/4))
          breaks <- make.even.breaks( diameter(window), bstep=eps/4)
        }

        return(breaks)
}
