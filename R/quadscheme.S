#
#
#      quadscheme.S
#
#      $Revision: 4.1 $    $Date: 2001/08/07 10:14:29 $
#
#      quadscheme()    generate a quadrature scheme from 
#		       data and dummy point patterns.
#
#      quadscheme.spatial()    case where both patterns are unmarked
#
#      quadscheme.replicated() case where data are multitype
#
#
#---------------------------------------------------------------------

quadscheme <- function(data, dummy=default.dummy(data), ...) {
        #
	# generate a quadrature scheme from data and dummy patterns.
	#
	# Other arguments control how the quadrature weights are computed
        #
  mX <- is.marked(data)
  mQ <- is.marked(dummy)
  
  if(!mX && !mQ)
    quadscheme.spatial(data, dummy, ...)
  else if(mX && !mQ)
    quadscheme.replicated(data, dummy, ...)
  else if(!mX && mQ)
    stop("dummy points are marked but data are unmarked")
  else
    stop("marked data and marked dummy points -- sorry, this case is not implemented")
}

quadscheme.spatial <-
  function(data, dummy=default.dummy(data), method="grid", ...) {
        #
	# generate a quadrature scheme from data and dummy patterns.
	#
	# The 'method' may be "grid" or "dirichlet"
	#
	# '...' are passed to gridweights() or dirichlet.weights()
        #
        # quadscheme.spatial:
        #       for unmarked point patterns.
        #
        #       weights are determined only by spatial locations
        #       (i.e. weight computations ignore any marks)
	#
        # No two points should have the same spatial location
        # 

	data <- as.ppp(data)
	dummy <- as.ppp(dummy, data$window)
		# note data$window is the DEFAULT quadrature window
		# unless otherwise specified in 'dummy'

        if(is.marked(data))
          warning("marks in data pattern - ignored")
        if(is.marked(dummy))
          warning("marks in dummy pattern - ignored")
        
	both <- as.ppp(concatxy(data, dummy), dummy$window)
	switch(method,
		grid={
			w <- gridweights(both, window= dummy$window, ...)
		},
		dirichlet = {
			w <- dirichlet.weights(both, window=dummy$window, ...)
		},
		{ 
			stop(paste("unrecognised method \'", method, "\'")) 
		}
	)
	Q <- quad(data, dummy, w)
	invisible(Q)
}

"quadscheme.replicated" <-
  function(data, dummy=default.dummy(data), method="grid", ...) {
        #
	# generate a quadrature scheme from data and dummy patterns.
	#
	# The 'method' may be "grid" or "dirichlet"
	#
	# '...' are passed to gridweights() or dirichlet.weights()
        #
        # quadscheme.replicated:
        #       for multitype point patterns.
        #
        # No two points in 'data'+'dummy' should have the same spatial location

	data <- as.ppp(data)
	dummy <- as.ppp(dummy, data$window)
		# note data$window is the DEFAULT quadrature window
		# unless otherwise specified in 'dummy'

        if(!is.marked(data))
          stop("data pattern does not have marks")
        if(is.marked(dummy))
          warning("dummy points have marks --- ignored")

        # first, ignore marks and compute spatial weights
        P <- quadscheme.spatial(unmark(data), dummy, method, ...)
        W <- w.quad(P)
        iz <- is.data(P)
        Wdat <- W[iz]
        Wdum <- W[!iz]

        # find the set of all possible marks

        if(!is.factor(data$marks))
          stop("data$marks is not a factor")
        markset <- levels(data$marks)
        nmarks <- length(markset)
        
        # replicate dummy points, one copy for each possible mark
        # -> dummy x {1,..,K}
        
        dumdum <- replicate(dummy, markset)
        Wdumdum <- rep(Wdum, nmarks)
        
        # also make dummy marked points at same locations as data points
        # but with different marks

        dumdat <- replicate(unmark(data), markset)
        Wdumdat <- rep(Wdat, nmarks)
        Mdumdat <- dumdat$marks
        
        Mrepdat <- rep(data$marks, nmarks)

        ok <- (Mdumdat != Mrepdat)
        dumdat <- dumdat[ok,]
        Wdumdat <- Wdumdat[ok]

        # combine the two dummy patterns

        dumb <- superimpose(dumdum, dumdat)
        Wdumb <- c(Wdumdum, Wdumdat)

        # wrap up

	Q <- quad(data, dumb, c(Wdat, Wdumb))
	invisible(Q)
}


"replicate" <-
function(pp, markset, fac=T) {
  # given an unmarked point pattern 'pp'
  # and a finite set of marks,
  # create the marked point pattern which is
  # the Cartesian product, consisting of all pairs (u,k)
  # where u is a point of 'pp' and k is a mark in 'markset'
  nmarks <- length(markset)
  result <- ppp(
                rep(pp$x, nmarks),
                rep(pp$y, nmarks),
                window=pp$window,
                marks=rep(markset, rep(pp$n, nmarks))
  )
  if(fac)
    result$marks <- factor(result$marks, levels=markset)
  result
}
