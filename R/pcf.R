#
#   pcf.R
#
#   $Revision: 1.6 $   $Date: 2004/01/13 12:12:35 $
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

        if(verifyclass(X, "fasp", fatal=FALSE)) {
          # function array - go to work on each function
          Y <- X
          Y$title <- paste("Array of pair correlation functions",
                           if(!is.null(X$dataname)) "for",
                           X$dataname)
          n <- length(X$fns)
          for(i in 1:n) {
            Xi <- X$fns[[i]]
            PCFi <- pcf(Xi, ..., method=method)
            Y$fns[[i]] <- as.fv(PCFi)
            if(is.fv(PCFi))
               Y$default.formula[[i]] <- attr(PCFi, "fmla")
          }
          return(Y)
        }

        if(is.fv(X)) {
          # extract r and the recommended estimate of K
          r <- X[[attr(X, "argu")]]
          K <- X[[attr(X, "valu")]]
          alim <- attr(X, "alim")
        } else if(inherits(X, "data.frame")) {
          # guess 
          r <- X$r
          K <- X$border
          alim <- NULL
        } else
          stop("X should be either a point pattern or the value returned by Kest() or Kcross() or alltypes(..., \"K\")")

	# remove NA's
	ok <- !is.na(K)
        K <- K[ok]
        r <- r[ok]
	switch(method,
		a = { 
			ss <- smooth.spline(r, K, ...)
			dK <- predict(ss, r, deriv=1)$y
			g <- dK/(2 * pi * r)
		},
		b = {
			y <- K/(2 * pi * r)
			y[is.nan(y)] <- 0
			ss <- smooth.spline(r, y, ...)
			dy <- predict(ss, r, deriv=1)$y
			g <- dy + y/r
		},
		c = {
			z <- K/(pi * r^2)
			z[is.nan(z)] <- 1
			ss <- smooth.spline(r, z, ...)
			dz <- predict(ss, r, deriv=1)$y
			g <- (r/2) * dz + z
		},
		stop(paste("unrecognised method \"", method, "\""))
	)

        # pack result into "fv" data frame
        Z <- fv(data.frame(r=r, pcf=g, theo=rep(1, length(r))),
                "r", "pcf(r)", "pcf", cbind(pcf, theo) ~ r, alim,
                c("r", "pcf(r)", "1"),
                c("distance argument r",
                  "estimate of pair correlation function pcf(r)",
                  "theoretical Poisson value, pcf(r) = 1"))
	return(Z)
}
