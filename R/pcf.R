#
#   pcf.R
#
#   $Revision: 1.2 $   $Date: 2001/11/23 06:25:02 $
#
#
#   calculate pair correlation function
#   from estimate of K or Kcross
#
#

"pcf" <-
function(X, ..., method="c") { 
	require(modreg)

	if(verifyclass(X, "ppp", fatal=FALSE))
        # point pattern - estimate K and continue
		X <- Kest(X)
        else if(verifyclass(X, "fasp", fatal=FALSE)) {
          # function array - go to work on each function
          n <- length(X$fns)
          for(i in 1:n) {
            Xi <- X$fns[[i]]
            if(!is.data.frame(Xi))
              stop("Internal error - an entry in the function array is not a data frame")
            X$fns[[i]] <- pcf(Xi, ..., method=method)
          }
	return(X)
        }
	else if(!is.data.frame(X) || is.null(X$r) || is.null(X$border))
		stop("X should be either a point pattern or the value returned by Kest() or Kcross() or alltypes(..., \"K\")")
        
	# remove NA's
	ok <- !is.na(X$border) 
	X <- X[ok, ]
	r <- X$r
	K <- X$border
	switch(method,
		a = { 
			ss <- smooth.spline(r, K, ...)
			dK <- predict.smooth.spline(ss, r, deriv=1)$y
			g <- dK/(2 * pi * r)
		},
		b = {
			y <- K/(2 * pi * r)
			y[is.nan(y)] <- 0
			ss <- smooth.spline(r, y, ...)
			dy <- predict.smooth.spline(ss, r, deriv=1)$y
			g <- dy + y/r
		},
		c = {
			z <- K/(pi * r^2)
			z[is.nan(z)] <- 1
			ss <- smooth.spline(r, z, ...)
			dz <- predict.smooth.spline(ss, r, deriv=1)$y
			g <- (r/2) * dz + z
		},
		stop(paste("unrecognised method \"", method, "\""))
	)
	X$pcf <- g
	return(X)
}
