#
#	$Revision: 1.5 $	$Date: 2005/04/29 20:28:48 $
#
#    ppm()
#          Fit a point process model to a two-dimensional point pattern
#
#

"ppm" <- 
function(Q,
         trend = ~1,
	 interaction = NULL,
         covariates = NULL,
	 correction="border",
	 rbord = 0,
         use.gam=FALSE,
         method = "mpl",
         forcefit=FALSE,
         nsim=100,
         nrmh=1e5,
         start=list(n.start=X$n),
         control=list(nrep=nrmh),
         verb=TRUE
) {
  if(!(method %in% c("mpl", "ho")))
    stop(paste("Unrecognised fitting method \"", method, "\"", sep=""))
  cl <- match.call()
  
  fitMPL <- mpl.engine(Q=Q, trend=trend,
                    interaction=interaction,
                    covariates=covariates,
                    correction=correction,
                    rbord=rbord, use.gam=use.gam,
                    forcefit=forcefit)
  fitMPL$call <- cl

  if(method == "mpl" || is.poisson.ppm(fitMPL))
    return(fitMPL)

  fitHO <- ho.engine(fitMPL, nsim=nsim, nrmh=nrmh, start=start,
                     control=control, verb=verb)

  if(!is.null(fitHO))
    return(fitHO)
  else
    return(fitMPL)

}

