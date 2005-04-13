#
#	$Revision: 1.3 $	$Date: 2005/04/12 20:23:57 $
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
         forcefit=FALSE
) {
  if(method != "mpl")
    stop(paste("Unrecognised fitting method \"", method, "\"", sep=""))

  fit <- mpl.engine(Q=Q, trend=trend,
                    interaction=interaction,
                    covariates=covariates,
                    correction=correction,
                    rbord=rbord, use.gam=use.gam,
                    forcefit=forcefit)
  fit$call <- deparse(sys.call())
  return(fit)
}

