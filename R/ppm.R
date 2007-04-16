#
#	$Revision: 1.12 $	$Date: 2007/03/28 02:59:44 $
#
#    ppm()
#          Fit a point process model to a two-dimensional point pattern
#
#

"ppm" <- 
function(Q,
         trend = ~1,
	 interaction = NULL,
         ..., 
         covariates = NULL,
	 correction="border",
	 rbord = 0,
         use.gam=FALSE,
         method = "mpl",
         forcefit=FALSE,
         nsim=100,
         nrmh=1e5,
         start=NULL,
         control=list(nrep=nrmh),
         verb=TRUE
) {
  if(!(method %in% c("mpl", "ho")))
    stop(paste("Unrecognised fitting method", sQuote(method)))
  cl <- match.call()
  callstring <- paste(deparse(sys.call()), collapse="")

  fitMPL <- mpl.engine(Q=Q, trend=trend,
                       interaction=interaction,
                       covariates=covariates,
                       correction=correction,
                       rbord=rbord, use.gam=use.gam,
                       forcefit=forcefit,
                       callstring=callstring,
                       preponly=FALSE,
                       ...)
  
  fitMPL$call <- cl
  fitMPL$callframe <- parent.frame()

  if(method == "mpl" || is.poisson.ppm(fitMPL))
    return(fitMPL)

  fitHO <- ho.engine(fitMPL, nsim=nsim, nrmh=nrmh, start=start,
                     control=control, verb=verb)

  if(!is.null(fitHO))
    return(fitHO)
  else
    return(fitMPL)

}

