#
#	$Revision: 1.2 $	$Date: 2004/06/09 10:58:18 $
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
         method = "mpl"
) {
  if(method != "mpl")
    stop(paste("Unrecognised fitting method \"", method, "\"", sep=""))

  fit <- mpl.engine(Q=Q, trend=trend,
                    interaction=interaction,
                    covariates=covariates,
                    correction=correction,
                    rbord=rbord, use.gam=use.gam)
  fit$call <- deparse(sys.call())
  return(fit)
}

